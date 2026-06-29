/*
 * Standalone decompressor for the LZ/Huffman format used by func_800157B8.
 *
 * Usage:
 *   decomp <input.bin> [output.bin]
 *
 * If output is omitted, writes to "<input>.out".
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define DECOMP_HUFFMAN_MAGIC   0x74
#define DECOMP_CTRL_END        0x00
#define DECOMP_CTRL_RLE_LONG   0x7F
#define DECOMP_CTRL_RLE_SHORT  0x7E
#define DECOMP_WINDOW_LIMIT    0x800
#define DECOMP_BYTE_MAP_SIZE   0x100
#define DECOMP_DEFAULT_OUT_CAP (4u * 1024u * 1024u)

typedef uint8_t  u8;
typedef uint16_t u16;
typedef int16_t  s16;
typedef uint32_t u32;
typedef int32_t  s32;

static u8* g_src;
static u8* g_src_start;
static u8* g_src_end;
static s32 g_pending_flag;
static u8  g_pending_byte;
static u8  g_byte_to_symbol[DECOMP_BYTE_MAP_SIZE];
static u8  g_symbol_code[256];   /* index 1..num_symbols -> raw code byte */
static u16 g_symbol_value[256];  /* index 1..num_symbols -> decoded u16 */
static u8  g_gap_dims_valid;
static u16 g_gap_w_units;
static u8  g_clut_uses_bytes;

static u8 Decomp_ReadHuffmanByte(void);
static s32 Decomp_Decompress(u8* src, size_t src_size, u8* dst, size_t dst_cap, size_t* out_size);
static s32 Decomp_DecompressAll(u8* src, size_t src_size, u8* dst, size_t dst_cap, size_t* out_size);

static u8 Decomp_ReadHuffmanByte(void) {
    u8* src;
    u8 map_index;
    u16 symbol;

    if (g_pending_flag != 0) {
        g_pending_flag = 0;
        return g_pending_byte;
    }

    if (g_src >= g_src_end) {
        fprintf(stderr, "error: unexpected end of input while reading Huffman stream\n");
        fprintf(stderr, "  offset: %td / %td bytes\n",
            (ptrdiff_t)(g_src - g_src_start), (ptrdiff_t)(g_src_end - g_src_start));
        exit(1);
    }

    src = g_src;
    map_index = g_byte_to_symbol[*src];

    if (map_index == 0) {
        g_src = src + 1;
        return *src;
    }

    symbol = g_symbol_value[map_index];
    g_pending_flag = 1;
    g_src = src + 1;
    g_pending_byte = (u8)(symbol >> 8);
    return (u8)symbol;
}

static u8* Decomp_CopyRepeat(u8* dst, u8* dst_end, u32 count, u8 value) {
    u32 i;

    for (i = 0; i < count; i++) {
        if (dst >= dst_end) {
            fprintf(stderr, "error: output buffer overflow during repeat copy\n");
            exit(1);
        }
        *dst++ = value;
    }
    return dst;
}

static u8* Decomp_CopyBackref(u8* dst, u8* dst_end, u8* base, u32 offset, u32 count) {
    u32 i;
    u8* ref_base;
    u8 use_safe_window;

    ref_base = dst;
    use_safe_window = (size_t)(ref_base - base) < DECOMP_WINDOW_LIMIT;

    for (i = 0; i < count; i++) {
        if (dst >= dst_end) {
            fprintf(stderr, "error: output buffer overflow during backref copy\n");
            exit(1);
        }
        if (use_safe_window && (size_t)(dst - ref_base) < offset) {
            *dst++ = 0;
        } else {
            *dst++ = *(dst - offset);
        }
    }
    return dst;
}

static u8* Decomp_CopyLiteralsHuffman(u8* dst, u8* dst_end, u32 count) {
    u32 i;

    for (i = 0; i < count; i++) {
        if (dst >= dst_end) {
            fprintf(stderr, "error: output buffer overflow during literal copy\n");
            exit(1);
        }
        *dst++ = Decomp_ReadHuffmanByte();
    }
    return dst;
}

static u8* Decomp_CopyLiteralsRaw(u8* dst, u8* dst_end, u8* src, u32 count) {
    u32 i;

    for (i = 0; i < count; i++) {
        if (dst >= dst_end) {
            fprintf(stderr, "error: output buffer overflow during raw literal copy\n");
            exit(1);
        }
        if (src >= g_src_end) {
            fprintf(stderr, "error: unexpected end of input during raw literal copy\n");
            exit(1);
        }
        *dst++ = *src++;
    }
    g_src = src;
    return dst;
}

static u32 Decomp_BackrefLength(u8 ctrl) {
    return ((ctrl >> 3) & 0xF) + 3;
}

static u32 Decomp_BackrefOffset(u8 ctrl, u8 offset_byte) {
    return offset_byte + ((ctrl << 8) & 0x700) + 1;
}

static u8* Decomp_HandleBlockHuffman(u8 ctrl, u8* dst, u8* dst_end, u8* base) {
    u32 count;
    u8 value;

    if (ctrl == DECOMP_CTRL_END) {
        return dst;
    }

    if (ctrl == DECOMP_CTRL_RLE_LONG) {
        value = Decomp_ReadHuffmanByte();
        count = value + ((u32)Decomp_ReadHuffmanByte() << 8);
        value = Decomp_ReadHuffmanByte();
        return Decomp_CopyRepeat(dst, dst_end, count, value);
    }

    if (ctrl == DECOMP_CTRL_RLE_SHORT) {
        count = Decomp_ReadHuffmanByte() + 4;
        value = Decomp_ReadHuffmanByte();
        return Decomp_CopyRepeat(dst, dst_end, count, value);
    }

    if (ctrl & 0x80) {
        count = Decomp_BackrefLength(ctrl);
        return Decomp_CopyBackref(dst, dst_end, base,
            Decomp_BackrefOffset(ctrl, Decomp_ReadHuffmanByte()), count);
    }

    return Decomp_CopyLiteralsHuffman(dst, dst_end, ctrl & 0x7F);
}

static void Decomp_InitHuffman(u8 num_symbols) {
    u32 byte_value;
    u32 symbol_index;

    memset(g_symbol_code, 0, sizeof(g_symbol_code));
    memset(g_symbol_value, 0, sizeof(g_symbol_value));
    memset(g_byte_to_symbol, 0, sizeof(g_byte_to_symbol));

    for (symbol_index = 1; symbol_index <= num_symbols; symbol_index++) {
        if (g_src >= g_src_end) {
            fprintf(stderr, "error: truncated Huffman symbol table\n");
            exit(1);
        }
        g_symbol_code[symbol_index] = *g_src++;
    }

    for (byte_value = 0; byte_value < DECOMP_BYTE_MAP_SIZE; byte_value++) {
        for (symbol_index = 1; symbol_index <= num_symbols; symbol_index++) {
            if (byte_value == g_symbol_code[symbol_index]) {
                g_byte_to_symbol[byte_value] = (u8)symbol_index;
                break;
            }
        }
    }

    for (symbol_index = 1; symbol_index <= num_symbols; symbol_index++) {
        if (g_src + 1 >= g_src_end) {
            fprintf(stderr, "error: truncated Huffman decode table\n");
            exit(1);
        }
        g_symbol_value[symbol_index] = (u16)(g_src[0] + (g_src[1] << 8));
        g_src += 2;
    }
}

static u8* Decomp_DecompressRaw(u8* dst, u8* dst_end, u8* base, u8* src_start) {
    u8* block;
    u8* end_ptr;
    u8 ctrl;
    u32 count;

    g_src = src_start + 1;
    end_ptr = src_start + 2;

    if (src_start + 1 >= g_src_end) {
        fprintf(stderr, "error: truncated raw header\n");
        exit(1);
    }

    if (src_start[1] == 0) {
        g_src = end_ptr;
        return dst;
    }

    do {
        if (g_src >= g_src_end) {
            fprintf(stderr, "error: unexpected end of input in raw block stream\n");
            exit(1);
        }

        block = g_src;
        ctrl = block[0];

        if (ctrl == DECOMP_CTRL_RLE_LONG) {
            if (block + 3 >= g_src_end) {
                fprintf(stderr, "error: truncated raw RLE long block\n");
                exit(1);
            }
            count = block[1] + ((u32)block[2] << 8);
            dst = Decomp_CopyRepeat(dst, dst_end, count, block[3]);
            g_src = block + 4;
        } else if (ctrl == DECOMP_CTRL_RLE_SHORT) {
            if (block + 2 >= g_src_end) {
                fprintf(stderr, "error: truncated raw RLE short block\n");
                exit(1);
            }
            count = block[1] + 4;
            dst = Decomp_CopyRepeat(dst, dst_end, count, block[2]);
            g_src = block + 3;
        } else if (ctrl & 0x80) {
            if (block + 1 >= g_src_end) {
                fprintf(stderr, "error: truncated raw backref block\n");
                exit(1);
            }
            count = Decomp_BackrefLength(ctrl);
            dst = Decomp_CopyBackref(dst, dst_end, base,
                Decomp_BackrefOffset(ctrl, block[1]), count);
            g_src = block + 2;
        } else {
            count = ctrl & 0x7F;
            if (block + count >= g_src_end) {
                fprintf(stderr, "error: truncated raw literal block\n");
                exit(1);
            }
            g_src = block + 1;
            dst = Decomp_CopyLiteralsRaw(dst, dst_end, g_src, count);
            g_src = block + 1 + count;
        }

        end_ptr = g_src + 1;
    } while (g_src < g_src_end && g_src[0] != 0);

    g_src = end_ptr;
    return dst;
}

static s32 Decomp_Decompress(u8* src, size_t src_size, u8* dst, size_t dst_cap, size_t* out_size) {
    u8* out;
    u8* base;
    u8* src_start;
    u8 ctrl;

    if (src_size == 0) {
        fprintf(stderr, "error: input file is empty\n");
        return -1;
    }

    src_start = src;
    g_src = src;
    g_src_start = src;
    g_src_end = src + src_size;
    g_pending_flag = 0;

    base = dst;
    out = dst;

    if (src[0] == DECOMP_HUFFMAN_MAGIC) {
        if (src_size < 2) {
            fprintf(stderr, "error: truncated Huffman header\n");
            return -1;
        }

        g_src = src + 2;
        Decomp_InitHuffman(src[1]);
        g_pending_flag = 0;

        for (;;) {
            ctrl = Decomp_ReadHuffmanByte();
            if (ctrl == DECOMP_CTRL_END) {
                break;
            }
            out = Decomp_HandleBlockHuffman(ctrl, out, dst + dst_cap, base);
        }
    } else {
        out = Decomp_DecompressRaw(out, dst + dst_cap, base, src);
    }

    *out_size = (size_t)(out - dst);
    return (s32)(g_src - src_start);
}

static int Decomp_IsGapHeader(const u8* src, size_t src_size, size_t src_pos) {
    if (src_pos + 4 > src_size) {
        return 0;
    }
    return src[src_pos] == 0x40 && src[src_pos + 1] == 0xFF &&
        src[src_pos + 2] == 0x00 && src[src_pos + 3] == 0x20;
}

static size_t Decomp_FindGapChunkOffset(const u8* src, size_t src_size, size_t src_pos) {
    size_t j;
    size_t limit;

    limit = src_pos + 128;
    if (limit > src_size) {
        limit = src_size;
    }

    for (j = src_pos; j + 1 < limit; j++) {
        if (src[j] == DECOMP_HUFFMAN_MAGIC) {
            return j;
        }
    }

    return src_size;
}

static void Decomp_ResetGapDims(void) {
    g_gap_dims_valid = 0;
    g_gap_w_units = 0;
}

static void Decomp_ParseGapHeader(const u8* src, size_t src_size, size_t src_pos) {
    u16 w_units;

    if (!Decomp_IsGapHeader(src, src_size, src_pos)) {
        return;
    }
    if (src_pos + 0x26 > src_size) {
        return;
    }

    w_units = (u16)(src[src_pos + 0x24] | (src[src_pos + 0x25] << 8));
    if (w_units > 0 && w_units <= 512) {
        g_gap_w_units = w_units;
        g_gap_dims_valid = 1;
    }
}

static int Decomp_TryResumeAfterGap(const u8* src, size_t src_size, size_t* src_pos, u8 used_gap) {
    size_t next_pos;

    if (used_gap || !Decomp_IsGapHeader(src, src_size, *src_pos)) {
        return 0;
    }

    Decomp_ParseGapHeader(src, src_size, *src_pos);
    next_pos = Decomp_FindGapChunkOffset(src, src_size, *src_pos);
    if (next_pos >= src_size || src[next_pos] != DECOMP_HUFFMAN_MAGIC) {
        return 0;
    }

    *src_pos = next_pos;
    return 1;
}

static s32 Decomp_DecompressAll(u8* src, size_t src_size, u8* dst, size_t dst_cap, size_t* out_size) {
    size_t src_pos;
    size_t total_out;
    u8* out;
    u8* out_end;
    u32 chunk_index;
    u8 used_gap;
    u8 gap_chunk_pending;

    if (src_size == 0) {
        fprintf(stderr, "error: input file is empty\n");
        return -1;
    }

    Decomp_ResetGapDims();
    src_pos = 0;
    total_out = 0;
    out = dst;
    out_end = dst + dst_cap;
    chunk_index = 0;
    used_gap = 0;
    gap_chunk_pending = 0;

    while (src_pos < src_size) {
        size_t chunk_out;
        s32 consumed;

        if (src[src_pos] == DECOMP_HUFFMAN_MAGIC) {
            consumed = Decomp_Decompress(src + src_pos, src_size - src_pos, out, (size_t)(out_end - out), &chunk_out);
        } else if (src_pos == 0) {
            consumed = Decomp_Decompress(src, src_size, out, (size_t)(out_end - out), &chunk_out);
        } else if (Decomp_TryResumeAfterGap(src, src_size, &src_pos, used_gap)) {
            used_gap = 1;
            gap_chunk_pending = 1;
            fprintf(stderr, "info: resuming after gap header, next chunk at offset %zu\n", src_pos);
            continue;
        } else {
            break;
        }

        if (consumed <= 0) {
            return -1;
        }

        chunk_index++;
        fprintf(stderr, "info: chunk %u at offset %zu, consumed %d, output %zu bytes\n",
            chunk_index, src_pos, consumed, chunk_out);

        out += chunk_out;
        total_out += chunk_out;
        src_pos += (size_t)consumed;

        if (gap_chunk_pending) {
            gap_chunk_pending = 0;
            break;
        }

        if (src_pos < src_size && src[src_pos] != DECOMP_HUFFMAN_MAGIC) {
            if (Decomp_TryResumeAfterGap(src, src_size, &src_pos, used_gap)) {
                used_gap = 1;
                gap_chunk_pending = 1;
                fprintf(stderr, "info: resuming after gap header, next chunk at offset %zu\n", src_pos);
                continue;
            }
            break;
        }
    }

    *out_size = total_out;
    return (s32)src_pos;
}

static u8* read_file(const char* path, size_t* size_out) {
    FILE* fp;
    long size;
    u8* data;

    fp = fopen(path, "rb");
    if (!fp) {
        perror(path);
        return NULL;
    }

    if (fseek(fp, 0, SEEK_END) != 0) {
        perror("fseek");
        fclose(fp);
        return NULL;
    }

    size = ftell(fp);
    if (size < 0) {
        perror("ftell");
        fclose(fp);
        return NULL;
    }

    if (fseek(fp, 0, SEEK_SET) != 0) {
        perror("fseek");
        fclose(fp);
        return NULL;
    }

    data = (u8*)malloc((size_t)size);
    if (!data) {
        fprintf(stderr, "error: out of memory\n");
        fclose(fp);
        return NULL;
    }

    if (fread(data, 1, (size_t)size, fp) != (size_t)size) {
        fprintf(stderr, "error: failed to read entire file\n");
        free(data);
        fclose(fp);
        return NULL;
    }

    fclose(fp);
    *size_out = (size_t)size;
    return data;
}

static int write_file(const char* path, const u8* data, size_t size) {
    FILE* fp;
    size_t written;

    fp = fopen(path, "wb");
    if (!fp) {
        perror(path);
        return 0;
    }

    written = fwrite(data, 1, size, fp);
    fclose(fp);

    if (written != size) {
        fprintf(stderr, "error: failed to write entire output file\n");
        return 0;
    }

    return 1;
}

static void make_default_output_name(const char* input, char* output, size_t output_cap) {
    const char* dot;
    size_t base_len;

    dot = strrchr(input, '.');
    if (dot != NULL) {
        base_len = (size_t)(dot - input);
        if (base_len + 5 < output_cap) {
            memcpy(output, input, base_len);
            memcpy(output + base_len, ".out", 5);
            return;
        }
    }

    if (strlen(input) + 5 < output_cap) {
        sprintf(output, "%s.out", input);
    } else {
        output[0] = '\0';
    }
}

static void psx_color_to_rgb(u16 c, u8* r, u8* g, u8* b) {
    *r = (u8)(((c >> 0) & 0x1F) << 3);
    *g = (u8)(((c >> 5) & 0x1F) << 3);
    *b = (u8)(((c >> 10) & 0x1F) << 3);
}

static int Decomp_ScorePixStart(const u8* data, size_t size, size_t pix_off) {
    size_t i;
    u8 seen[256];
    int score;
    size_t unique;

    if (pix_off >= size) {
        return -999;
    }

    memset(seen, 0, sizeof(seen));
    unique = 0;
    for (i = 0; i < 64 && pix_off + i < size; i++) {
        if (!seen[data[pix_off + i]]) {
            seen[data[pix_off + i]] = 1;
            unique++;
        }
    }

    score = (int)unique;

    if (pix_off + 16 <= size) {
        int all_84 = 1;
        for (i = 0; i < 16; i++) {
            if (data[pix_off + i] != 0x84) {
                all_84 = 0;
                break;
            }
        }
        if (all_84) {
            score -= 100;
        }
    }

    if (pix_off + 32 <= size) {
        int all_zero = 1;
        int has_later = 0;
        size_t i;

        for (i = 0; i < 32; i++) {
            if (data[pix_off + i] != 0) {
                all_zero = 0;
                break;
            }
        }
        if (all_zero) {
            for (i = 0; i < 256 && pix_off + i < size; i++) {
                if (data[pix_off + i] != 0) {
                    has_later = 1;
                    break;
                }
            }
            score -= has_later ? 5 : 20;
        }
    }

    return score;
}

static size_t Decomp_ImgBlockOffset(const u8* data, size_t size) {
    u32 flags;
    u32 block_len;
    size_t img_words;
    size_t img_bytes;
    size_t pix_words;
    size_t pix_bytes;
    int score_words;
    int score_bytes;

    if (size < 20 || data[0] != 0x10) {
        return 0;
    }

    flags = (u32)data[4] | ((u32)data[5] << 8) | ((u32)data[6] << 16) | ((u32)data[7] << 24);
    if ((flags & 8) == 0) {
        return 8;
    }

    block_len = (u32)data[8] | ((u32)data[9] << 8) | ((u32)data[10] << 16) | ((u32)data[11] << 24);

    img_words = 8 + (size_t)block_len * 4;
    img_bytes = 8 + (size_t)block_len;
    pix_words = img_words + 12;
    pix_bytes = img_bytes + 12;

    score_words = pix_words < size ? Decomp_ScorePixStart(data, size, pix_words) : -999;
    score_bytes = pix_bytes < size ? Decomp_ScorePixStart(data, size, pix_bytes) : -999;

    g_clut_uses_bytes = 0;
    if (block_len == 44 || block_len == 524) {
        g_clut_uses_bytes = 1;
        return img_bytes;
    }
    if (score_bytes > score_words) {
        g_clut_uses_bytes = 1;
        return img_bytes;
    }
    return img_words;
}

static int Decomp_ClutBlockUsesBytes(const u8* data, size_t img_off) {
    (void)data;
    (void)img_off;
    return g_clut_uses_bytes;
}

static void Decomp_StandardizeClutBlockLen(u8* tim_buf, size_t img_off) {
    u32 block_len;
    u32 block_bytes;

    if (g_clut_uses_bytes || img_off < 12) {
        return;
    }

    block_len = (u32)tim_buf[8] | ((u32)tim_buf[9] << 8) | ((u32)tim_buf[10] << 16) | ((u32)tim_buf[11] << 24);
    if (block_len == 44 || block_len == 524) {
        return;
    }

    block_bytes = block_len * 4;
    tim_buf[8] = (u8)(block_bytes);
    tim_buf[9] = (u8)(block_bytes >> 8);
    tim_buf[10] = (u8)(block_bytes >> 16);
    tim_buf[11] = (u8)(block_bytes >> 24);
}

static int Decomp_ScoreDims(u32 bpp_mode, u16 w_units, u16 height);

static int Decomp_TryFitRawDims(u32 bpp_mode, size_t raw_pix, u16* w_units, u16* height,
    size_t* pix_bytes, size_t* trim_tail) {
    static const u16 width_units[] = { 160, 128, 80, 64, 40, 256, 96, 48, 32, 16, 320, 512 };
    size_t i;
    u16 wu;
    u16 h;
    size_t row_bytes;
    size_t expected;
    int best_score;
    u16 best_wu;
    u16 best_h;
    size_t best_expected;
    size_t best_trim;

    if (bpp_mode != 0 && bpp_mode != 1) {
        return 0;
    }

    best_score = -1;
    best_wu = 0;
    best_h = 0;
    best_expected = 0;
    best_trim = 0;

    for (i = 0; i < sizeof(width_units) / sizeof(width_units[0]); i++) {
        int score;
        size_t trim;

        wu = width_units[i];
        row_bytes = (size_t)wu * 2;
        if (row_bytes == 0 || raw_pix < row_bytes) {
            continue;
        }

        h = (u16)(raw_pix / row_bytes);
        if (h == 0) {
            continue;
        }

        expected = (size_t)wu * h * 2;
        if (expected == raw_pix) {
            score = Decomp_ScoreDims(bpp_mode, wu, h) + 20;
            trim = 0;
        } else if (expected < raw_pix && raw_pix - expected < row_bytes) {
            score = Decomp_ScoreDims(bpp_mode, wu, h) + 10;
            trim = raw_pix - expected;
        } else {
            continue;
        }

        if (score > best_score) {
            best_score = score;
            best_wu = wu;
            best_h = h;
            best_expected = expected;
            best_trim = trim;
        }
    }

    if (best_score < 0) {
        return 0;
    }

    *w_units = best_wu;
    *height = best_h;
    *pix_bytes = best_expected;
    *trim_tail = best_trim;
    return 1;
}

static u16 Decomp_WidthPxFromUnits(u32 bpp_mode, u16 w_units) {
    if (bpp_mode == 0) {
        return (u16)(w_units * 4);
    }
    if (bpp_mode == 1) {
        return (u16)(w_units * 2);
    }
    return w_units;
}

static size_t Decomp_ExpectedPixBytes(u16 w_units, u16 height) {
    return (size_t)w_units * height * 2;
}

static int Decomp_ValidateTimDims(u16 w_units, u16 height, size_t raw_pix, size_t pix_off) {
    size_t expected;

    if (w_units == 0 || height == 0) {
        return 0;
    }

    expected = Decomp_ExpectedPixBytes(w_units, height);
    if (expected == raw_pix) {
        return 1;
    }
    if (expected == raw_pix + pix_off) {
        return 2;
    }
    return 0;
}

static int Decomp_ParseFooterDims(const u8* footer, size_t footer_len, u16* w_units, u16* height,
    size_t raw_pix, size_t pix_off) {
    u16 wu;
    u16 h;
    size_t off;

    static const size_t footer_offsets[][2] = {
        { 0x0C, 0x12 },
        { 0x10, 0x12 },
        { 0x0C, 0x10 },
    };
    size_t i;

    if (footer == NULL || footer_len < 0x14) {
        return 0;
    }

    for (i = 0; i < sizeof(footer_offsets) / sizeof(footer_offsets[0]); i++) {
        off = footer_offsets[i][0];
        if (off + 2 > footer_len || footer_offsets[i][1] + 2 > footer_len) {
            continue;
        }
        wu = (u16)(footer[off] | (footer[off + 1] << 8));
        off = footer_offsets[i][1];
        h = (u16)(footer[off] | (footer[off + 1] << 8));
        if (Decomp_ValidateTimDims(wu, h, raw_pix, pix_off)) {
            *w_units = wu;
            *height = h;
            return 1;
        }
    }

    return 0;
}

static int Decomp_ParseStoredImageDims(const u8* data, size_t img_off, u16* w_units, u16* height,
    size_t raw_pix, size_t pix_off) {
    u16 wu;
    u16 h;

    wu = (u16)(data[img_off + 8] | (data[img_off + 9] << 8));
    h = (u16)(data[img_off + 10] | (data[img_off + 11] << 8));
    if (Decomp_ValidateTimDims(wu, h, raw_pix, pix_off)) {
        *w_units = wu;
        *height = h;
        return 1;
    }
    return 0;
}

static int Decomp_ScoreDims(u32 bpp_mode, u16 w_units, u16 height) {
    u16 width_px;
    int score;

    width_px = Decomp_WidthPxFromUnits(bpp_mode, w_units);
    if (width_px < 16 || width_px > 1024 || height < 1 || height > 1024) {
        return -1;
    }

    score = 0;
    if ((width_px % 16) == 0) {
        score += 4;
    }
    if (width_px == 256 || width_px == 128 || width_px == 320 || width_px == 512) {
        score += 3;
    }
    if (height == 256 || height == 128 || height == 64) {
        score += 2;
    }
    if ((width_px & (width_px - 1)) == 0) {
        score += 1;
    }
    if (height > 0) {
        u32 aspect = width_px > height ? width_px / height : height / width_px;
        if (aspect <= 4) {
            score += 2;
        }
    }
    if (width_px < 128) {
        score -= 15;
    }
    if (height > 256) {
        score -= 15;
    }
    return score;
}

static int Decomp_ScanFooterDims(const u8* footer, size_t footer_len, u32 bpp_mode,
    u16* w_units, u16* height, size_t raw_pix, size_t pix_off) {
    size_t w_off;
    size_t h_off;
    int best_score;
    u16 best_wu;
    u16 best_h;

    if (footer == NULL || footer_len < 4) {
        return 0;
    }

    best_score = -1;
    best_wu = 0;
    best_h = 0;

    for (w_off = 0; w_off + 2 <= footer_len; w_off += 2) {
        u16 wu = (u16)(footer[w_off] | (footer[w_off + 1] << 8));
        if (wu == 0 || wu > 512) {
            continue;
        }
        for (h_off = 0; h_off + 2 <= footer_len; h_off += 2) {
            u16 h = (u16)(footer[h_off] | (footer[h_off + 1] << 8));
            int score;

            if (w_off == h_off || h == 0 || h > 1024) {
                continue;
            }
            if (!Decomp_ValidateTimDims(wu, h, raw_pix, pix_off)) {
                continue;
            }

            score = Decomp_ScoreDims(bpp_mode, wu, h);
            if (w_off == 0x0C && h_off == 0x12) {
                score += 10;
            }
            if (score > best_score) {
                best_score = score;
                best_wu = wu;
                best_h = h;
            }
        }
    }

    if (best_score < 0) {
        return 0;
    }

    *w_units = best_wu;
    *height = best_h;
    return 1;
}

static int Decomp_InferDimsFromSize(u32 bpp_mode, size_t raw_pix, size_t pix_off,
    u16* w_units, u16* height) {
    size_t totals[2];
    size_t total_index;
    int best_score;
    u16 best_wu;
    u16 best_h;
    u16 wu;
    u16 h;
    size_t product;

    totals[0] = raw_pix;
    totals[1] = raw_pix + pix_off;
    best_score = -1;
    best_wu = 0;
    best_h = 0;

    for (total_index = 0; total_index < 2; total_index++) {
        size_t total;

        total = totals[total_index];
        if ((total % 2) != 0) {
            continue;
        }

        product = total / 2;
        for (wu = 1; wu <= 512; wu++) {
            if ((product % wu) != 0) {
                continue;
            }
            h = (u16)(product / wu);
            if (h == 0) {
                continue;
            }
            if (!Decomp_ValidateTimDims(wu, h, raw_pix, pix_off)) {
                continue;
            }
            {
                int score = Decomp_ScoreDims(bpp_mode, wu, h);
                if (score > best_score) {
                    best_score = score;
                    best_wu = wu;
                    best_h = h;
                }
            }
        }
    }

    if (best_score < 0) {
        return 0;
    }

    *w_units = best_wu;
    *height = best_h;
    return 1;
}

static int Decomp_ParseGapOutputDims(size_t size, size_t raw_pix, size_t pix_off,
    u16* w_units, u16* height) {
    u16 h;

    if (!g_gap_dims_valid || g_gap_w_units == 0) {
        return 0;
    }

    if ((size % ((size_t)g_gap_w_units * 2)) != 0) {
        return 0;
    }

    h = (u16)(size / ((size_t)g_gap_w_units * 2));
    if (h == 0) {
        return 0;
    }

    if (!Decomp_ValidateTimDims(g_gap_w_units, h, raw_pix, pix_off)) {
        return 0;
    }

    *w_units = g_gap_w_units;
    *height = h;
    return 1;
}

static int Decomp_ResolveTimDims(const u8* data, size_t size, const u8* footer, size_t footer_len,
    u32 bpp_mode, u16* w_units, u16* height, size_t* pix_bytes, size_t* lead_pad, size_t* trim_tail) {
    size_t clut_end;
    size_t img_off;
    size_t pix_off;
    size_t raw_pix;
    size_t expected_pix;
    int clut_bytes;

    clut_end = Decomp_ImgBlockOffset(data, size);
    if (clut_end == 0 || clut_end + 12 > size) {
        return 0;
    }

    img_off = clut_end;
    pix_off = img_off + 12;
    if (pix_off > size) {
        return 0;
    }

    raw_pix = size - pix_off;
    if (raw_pix == 0) {
        return 0;
    }

    *lead_pad = 0;
    *trim_tail = 0;
    *pix_bytes = 0;
    clut_bytes = Decomp_ClutBlockUsesBytes(data, img_off);
    {
        int dims_resolved = 0;

        if (clut_bytes && bpp_mode == 1) {
            u32 block_len = (u32)data[8] | ((u32)data[9] << 8) | ((u32)data[10] << 16) | ((u32)data[11] << 24);
            if (block_len == 524) {
                u16 h = (u16)(raw_pix / 320);
                if (h > 0 && (size_t)h * 320 <= raw_pix && raw_pix - (size_t)h * 320 < 320) {
                    *w_units = 160;
                    *height = h;
                    *pix_bytes = raw_pix;
                    *trim_tail = 0;
                    dims_resolved = 1;
                }
            }
        }

        if (!dims_resolved && clut_bytes &&
            Decomp_TryFitRawDims(bpp_mode, raw_pix, w_units, height, pix_bytes, trim_tail)) {
            dims_resolved = 1;
        } else if (!dims_resolved && Decomp_ParseGapOutputDims(size, raw_pix, pix_off, w_units, height)) {
            dims_resolved = 1;
        } else if (!dims_resolved &&
            Decomp_TryFitRawDims(bpp_mode, raw_pix, w_units, height, pix_bytes, trim_tail)) {
            dims_resolved = 1;
        } else if (!dims_resolved && Decomp_ParseFooterDims(footer, footer_len, w_units, height, raw_pix, pix_off)) {
            dims_resolved = 1;
        } else if (!dims_resolved &&
            Decomp_ScanFooterDims(footer, footer_len, bpp_mode, w_units, height, raw_pix, pix_off)) {
            dims_resolved = 1;
        } else if (!dims_resolved &&
            Decomp_ParseStoredImageDims(data, img_off, w_units, height, raw_pix, pix_off)) {
            dims_resolved = 1;
        } else if (!dims_resolved && !Decomp_InferDimsFromSize(bpp_mode, raw_pix, pix_off, w_units, height)) {
            fprintf(stderr, "warning: export skipped - could not determine texture size "
                "(raw pixels %zu, header gap %zu)\n", raw_pix, pix_off);
            return 0;
        }
    }

    if (*pix_bytes == 0) {
        expected_pix = Decomp_ExpectedPixBytes(*w_units, *height);
        if (expected_pix == raw_pix) {
            *pix_bytes = expected_pix;
        } else if (!clut_bytes && expected_pix == raw_pix + pix_off) {
            *lead_pad = pix_off;
            *pix_bytes = expected_pix;
        } else {
            fprintf(stderr, "warning: export skipped - size mismatch for %u x %u units "
                "(expected %zu bytes, have %zu + %zu pad)\n",
                *w_units, *height, expected_pix, raw_pix, pix_off);
            return 0;
        }
    }

    return 1;
}

static void Decomp_WriteU32Le(FILE* fp, u32 value) {
    fputc((int)(value & 0xFF), fp);
    fputc((int)((value >> 8) & 0xFF), fp);
    fputc((int)((value >> 16) & 0xFF), fp);
    fputc((int)((value >> 24) & 0xFF), fp);
}

static int Decomp_ExportTimAndBmp(const u8* data, size_t size, const char* base_path,
    const u8* footer, size_t footer_len) {
    char tim_path[4096];
    char bmp_path[4096];
    size_t clut_end;
    size_t img_off;
    size_t pix_off;
    size_t pix_bytes;
    size_t lead_pad;
    size_t trim_tail;
    size_t tim_size;
    u32 img_words;
    u32 bpp_mode;
    u32 flags;
    u16 w_units;
    u16 height;
    u16 width;
    u8* tim_buf;
    FILE* fp;
    size_t x;
    size_t y;
    size_t row_size;
    size_t row_pad;
    u32 bmp_size;
    u8 palette[256 * 4];
    const u16* clut;
    size_t clut_count;
    size_t clut_words;
    size_t raw_pix;

    if (size < 16) {
        fprintf(stderr, "warning: export skipped - output too small (%zu bytes)\n", size);
        return 0;
    }

    if (data[0] != 0x10) {
        fprintf(stderr, "warning: export skipped - no TIM magic (first byte 0x%02X)\n", data[0]);
        return 0;
    }

    flags = (u32)data[4] | ((u32)data[5] << 8) | ((u32)data[6] << 16) | ((u32)data[7] << 24);
    bpp_mode = flags & 7;
    if (bpp_mode != 0 && bpp_mode != 1) {
        fprintf(stderr, "warning: export skipped - unsupported TIM bpp mode %u\n", bpp_mode);
        return 0;
    }

    clut_end = Decomp_ImgBlockOffset(data, size);
    if (clut_end == 0 || clut_end + 12 > size) {
        fprintf(stderr, "warning: export skipped - invalid CLUT/image layout\n");
        return 0;
    }

    img_off = clut_end;
    pix_off = img_off + 12;
    raw_pix = size - pix_off;

    if (!Decomp_ResolveTimDims(data, size, footer, footer_len, bpp_mode,
            &w_units, &height, &pix_bytes, &lead_pad, &trim_tail)) {
        return 0;
    }

    width = Decomp_WidthPxFromUnits(bpp_mode, w_units);
    img_words = (u32)(12 + pix_bytes);
    tim_size = img_off + 12 + pix_bytes;
    tim_buf = (u8*)malloc(tim_size);
    if (!tim_buf) {
        fprintf(stderr, "warning: export skipped - out of memory\n");
        return 0;
    }

    memcpy(tim_buf, data, img_off);
    Decomp_StandardizeClutBlockLen(tim_buf, img_off);
    tim_buf[img_off + 0] = (u8)(img_words);
    tim_buf[img_off + 1] = (u8)(img_words >> 8);
    tim_buf[img_off + 2] = (u8)(img_words >> 16);
    tim_buf[img_off + 3] = (u8)(img_words >> 24);
    tim_buf[img_off + 4] = data[img_off + 4];
    tim_buf[img_off + 5] = data[img_off + 5];
    tim_buf[img_off + 6] = data[img_off + 6];
    tim_buf[img_off + 7] = data[img_off + 7];
    tim_buf[img_off + 8] = (u8)(w_units);
    tim_buf[img_off + 9] = (u8)(w_units >> 8);
    tim_buf[img_off + 10] = (u8)(height);
    tim_buf[img_off + 11] = (u8)(height >> 8);

    if (lead_pad > 0) {
        memset(tim_buf + pix_off, 0, lead_pad);
        memcpy(tim_buf + pix_off + lead_pad, data + pix_off, pix_bytes);
    } else if (trim_tail > 0) {
        memcpy(tim_buf + pix_off, data + pix_off, pix_bytes);
    } else {
        memcpy(tim_buf + pix_off, data + pix_off, pix_bytes);
    }

    if (snprintf(tim_path, sizeof(tim_path), "%s.tim", base_path) >= (int)sizeof(tim_path)) {
        free(tim_buf);
        return 0;
    }
    if (!write_file(tim_path, tim_buf, tim_size)) {
        fprintf(stderr, "warning: failed to write %s\n", tim_path);
        free(tim_buf);
        return 0;
    }

    clut = (const u16*)(data + 20);
    if (g_clut_uses_bytes) {
        u32 block_len = (u32)data[8] | ((u32)data[9] << 8) | ((u32)data[10] << 16) | ((u32)data[11] << 24);
        clut_count = block_len >= 12 ? ((size_t)block_len - 12) / 2 : 0;
    } else {
        clut_words = (u32)data[8] | ((u32)data[9] << 8) | ((u32)data[10] << 16) | ((u32)data[11] << 24);
        clut_count = ((size_t)clut_words * 4 - 12) / 2;
    }
    if (clut_count > 256) {
        clut_count = 256;
    }

    for (x = 0; x < 256; x++) {
        u8 r;
        u8 g;
        u8 b;
        u16 color;

        if (x < clut_count) {
            color = clut[x];
        } else {
            color = 0;
        }
        psx_color_to_rgb(color, &r, &g, &b);
        palette[x * 4 + 0] = b;
        palette[x * 4 + 1] = g;
        palette[x * 4 + 2] = r;
        palette[x * 4 + 3] = 0;
    }

    row_size = (bpp_mode == 0) ? (size_t)(width / 2) : (size_t)width;
    row_pad = (4 - (row_size % 4)) % 4;
    bmp_size = 14 + 40 + 256 * 4 + (u32)((row_size + row_pad) * height);

    if (snprintf(bmp_path, sizeof(bmp_path), "%s.bmp", base_path) >= (int)sizeof(bmp_path)) {
        free(tim_buf);
        return 0;
    }

    fp = fopen(bmp_path, "wb");
    if (!fp) {
        fprintf(stderr, "warning: failed to write %s\n", bmp_path);
        free(tim_buf);
        return 0;
    }

    fputc('B', fp);
    fputc('M', fp);
    Decomp_WriteU32Le(fp, bmp_size);
    Decomp_WriteU32Le(fp, 0);
    Decomp_WriteU32Le(fp, 54 + 256 * 4);
    Decomp_WriteU32Le(fp, 40);
    Decomp_WriteU32Le(fp, (u32)width);
    Decomp_WriteU32Le(fp, (u32)height);
    fputc(1, fp);
    fputc(0, fp);
    fputc(8, fp);
    fputc(0, fp);
    Decomp_WriteU32Le(fp, (u32)((row_size + row_pad) * height));
    Decomp_WriteU32Le(fp, 0);
    Decomp_WriteU32Le(fp, 0);
    Decomp_WriteU32Le(fp, 256);
    Decomp_WriteU32Le(fp, 0);
    fwrite(palette, 1, sizeof(palette), fp);

    for (y = 0; y < height; y++) {
        const u8* row = tim_buf + pix_off + (height - 1 - y) * row_size;
        if (bpp_mode == 0) {
            fwrite(row, 1, row_size, fp);
        } else {
            for (x = 0; x < width; x++) {
                fputc(row[x], fp);
            }
        }
        for (x = 0; x < row_pad; x++) {
            fputc(0, fp);
        }
    }

    fclose(fp);
    free(tim_buf);

    printf("  exported TIM: %s (%u x %u, %s)\n", tim_path, width, height,
        bpp_mode == 0 ? "4bpp" : "8bpp");
    printf("  exported BMP: %s\n", bmp_path);
    return 1;
}

static void Decomp_MakeBaseName(const char* path, char* base, size_t base_cap) {
    const char* name;
    const char* dot;
    size_t len;

    name = strrchr(path, '\\');
    if (name == NULL) {
        name = strrchr(path, '/');
    }
    if (name == NULL) {
        name = path;
    } else {
        name++;
    }

    dot = strrchr(name, '.');
    if (dot != NULL) {
        len = (size_t)(dot - name);
    } else {
        len = strlen(name);
    }

    if (len + 1 > base_cap) {
        base[0] = '\0';
        return;
    }

    memcpy(base, name, len);
    base[len] = '\0';
}

static void print_usage(const char* prog) {
    fprintf(stderr, "Usage: %s <input.bin> [output.bin]\n", prog);
    fprintf(stderr, "\nDecompresses PS1 LZ/Huffman archives (func_800157B8 format).\n");
    fprintf(stderr, "Also writes <output>.tim and <output>.bmp when the data is a TIM texture.\n");
    fprintf(stderr, "If output is omitted, writes to \"<input>.out\".\n");
}

int main(int argc, char** argv) {
    const char* input_path;
    char default_output[4096];
    char input_base[4096];
    const char* output_path;
    u8* input_data;
    u8* output_data;
    size_t input_size;
    size_t output_size;
    size_t output_cap;
    s32 consumed;

    if (argc < 2 || argc > 3) {
        print_usage(argv[0]);
        return 1;
    }

    input_path = argv[1];
    if (argc == 3) {
        output_path = argv[2];
    } else {
        make_default_output_name(input_path, default_output, sizeof(default_output));
        if (default_output[0] == '\0') {
            fprintf(stderr, "error: output path too long\n");
            return 1;
        }
        output_path = default_output;
    }

    input_data = read_file(input_path, &input_size);
    if (!input_data) {
        return 1;
    }

    output_cap = input_size * 16;
    if (output_cap < DECOMP_DEFAULT_OUT_CAP) {
        output_cap = DECOMP_DEFAULT_OUT_CAP;
    }

    output_data = (u8*)malloc(output_cap);
    if (!output_data) {
        fprintf(stderr, "error: out of memory\n");
        free(input_data);
        return 1;
    }

    consumed = Decomp_DecompressAll(input_data, input_size, output_data, output_cap, &output_size);
    if (consumed < 0) {
        free(output_data);
        free(input_data);
        return 1;
    }

    if (!write_file(output_path, output_data, output_size)) {
        free(output_data);
        free(input_data);
        return 1;
    }

    if (!Decomp_ExportTimAndBmp(output_data, output_size, output_path,
            input_data + (size_t)consumed, input_size - (size_t)consumed)) {
        fprintf(stderr, "hint: check footer metadata or TIM header for texture dimensions.\n");
    }

    Decomp_MakeBaseName(input_path, input_base, sizeof(input_base));
    if (input_base[0] != '\0') {
        char alt_path[4096];
        const char* slash;
        size_t dir_len;

        slash = strrchr(output_path, '\\');
        if (slash == NULL) {
            slash = strrchr(output_path, '/');
        }

        if (slash != NULL) {
            dir_len = (size_t)(slash + 1 - output_path);
            if (snprintf(alt_path, sizeof(alt_path), "%.*s%s", (int)dir_len, output_path, input_base)
                < (int)sizeof(alt_path)) {
                if (strcmp(alt_path, output_path) != 0) {
                    Decomp_ExportTimAndBmp(output_data, output_size, alt_path,
                        input_data + (size_t)consumed, input_size - (size_t)consumed);
                }
            }
        } else if (strcmp(input_base, output_path) != 0) {
            Decomp_ExportTimAndBmp(output_data, output_size, input_base,
                input_data + (size_t)consumed, input_size - (size_t)consumed);
        }
    }

    printf("Decompressed %s -> %s\n", input_path, output_path);
    printf("  input size:  %zu bytes\n", input_size);
    printf("  output size: %zu bytes\n", output_size);
    printf("  compressed used: %d bytes\n", consumed);
    if ((size_t)consumed < input_size) {
        printf("  trailing data:   %zu bytes (file footer/metadata)\n", input_size - (size_t)consumed);
    }

    free(output_data);
    free(input_data);
    return 0;
}
