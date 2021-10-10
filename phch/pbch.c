/*
 * pbch.c
 *
 *  Created on: Jun 25, 2020
 *      Author: administrator
 */

#include "srslte/srslte.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "srslte/phy/common/phy_common.h"
#include "srslte/phy/phch/pbch.h"
#include "srslte/phy/utils/bit.h"
#include "srslte/phy/utils/debug.h"
#include "srslte/phy/utils/vector.h"
#include "srslte/phy/ue/ue_config.h"

bool srslte_pbch_exists(int nframe, int nslot, srslte_ssb_type_t ssb_type) {
	return false;
}

int srslte_pbch_res_init(srslte_grid_res_t* q, srslte_cell_t* cell) {
	uint16_t offset[SRSLTE_MAX_PORTS] = { 0 };
	srslte_grid_symbol_t symbol[SRSLTE_N_SLOT_SYMB] = { NULL_SYMB };

	offset[0] = cell->id % 4;
	srslte_grid_res_init(q, offset, 3, 1, 1, cell->sz);

	symbol[1] = RS_SYMB;
	srslte_grid_add_block(q, symbol, cell->ssb.ssb_offset, SRSLTE_PBCH_RB_NUM, 0);

	symbol[1] = NULL_SYMB;
	symbol[2] = RS_SYMB;
	srslte_grid_add_block(q, symbol, cell->ssb.ssb_offset, 4, 1);
	srslte_grid_add_block(q, symbol, cell->ssb.ssb_offset + 192, 4, 1);

	symbol[2] = NULL_SYMB;
	symbol[3] = RS_SYMB;
	srslte_grid_add_block(q, symbol, cell->ssb.ssb_offset, SRSLTE_PBCH_RB_NUM, 1);
	return SRSLTE_SUCCESS;
}

int srslte_pbch_put(cf_t* pbch, cf_t* slot_data, srslte_grid_res_t* q) {
	return srslte_grid_cp_re(q, pbch, slot_data, 0, true);
}

int srslte_pbch_get(cf_t* slot_data, cf_t* pbch, srslte_grid_res_t* q) {
	return srslte_grid_cp_re(q, slot_data, pbch, 0, false);
}

int srslte_pbch_ce_put(cf_t* ce, cf_t* ce_slot, srslte_grid_res_t* q) {
	return srslte_grid_cp_ce(q, ce, ce_slot, 0, true);
}

int srslte_pbch_ce_get(cf_t* ce_slot, cf_t* ce, srslte_grid_res_t* q) {
	return srslte_grid_cp_ce(q, ce_slot, ce, 0, false);
}

int srslte_pbch_refs_put(cf_t* pbch, cf_t* slot_data, srslte_grid_res_t* q) {
	return srslte_grid_cp_rs(q, pbch, slot_data, 0, true);
}

int srslte_pbch_refs_get(cf_t* slot_data, cf_t* pbch, srslte_grid_res_t* q) {
	return srslte_grid_cp_rs(q, slot_data, pbch, 0, false);
}

int srslte_pbch_dmrs_seq(srslte_sequence_t* seq, uint32_t cell_id, uint8_t ssb_index, uint8_t Lmax, uint8_t n_hf) {
	uint32_t iprime_ssb, c_init;

	if (Lmax == 4) {
		iprime_ssb = (ssb_index & 0x03) + (n_hf << 2);
	} else {
		iprime_ssb = ssb_index & 0x07;
	}
	c_init = ((iprime_ssb + 1) * (cell_id / 4 + 1) << 11) + ((iprime_ssb + 1) << 6) + (cell_id % 4);

	return srslte_sequence_LTE_pr(seq, 2 * SRSLTE_PBCH_DMRS_LEN, c_init);
}

int srslte_pbch_dmrs_gen(cf_t* dmrs, srslte_sequence_t* seq, uint32_t cell_id, uint8_t ssb_index, uint8_t Lmax,
		uint8_t n_hf) {
	srslte_pbch_dmrs_seq(seq, cell_id, ssb_index, Lmax, n_hf);
	for (int i = 0; i < SRSLTE_PBCH_DMRS_LEN; i++) {
		__real__ dmrs[i] = M_SQRT1_2 * (1 - 2 * seq->c[2 * i]);
		__imag__ dmrs[i] = M_SQRT1_2 * (1 - 2 * seq->c[2 * i + 1]);
	}
	return SRSLTE_PBCH_DMRS_LEN;
}

float srslte_pbch_refs_corr(srslte_chest_dl_t* q, cf_t* input[SRSLTE_MAX_PORTS]) {
	int n_hf, i_ssb;
	int hf_pos, ssb_pos;
	srslte_sequence_t seq;
	cf_t peak_value = 0, peak_max = 0;
	cf_t pilot_recv_signal[SRSLTE_PBCH_DMRS_LEN];
	cf_t pilot_refs_signal[2][8][SRSLTE_PBCH_DMRS_LEN];

	bzero(&seq, sizeof(srslte_sequence_t));
	srslte_pbch_refs_get(input[0], pilot_recv_signal, &q->res);

	for (n_hf = 0; n_hf < 2; n_hf++) {
		for (i_ssb = 0; i_ssb < 8; i_ssb++) {
			srslte_pbch_dmrs_gen(pilot_refs_signal[n_hf][i_ssb], &seq, q->cell.id, i_ssb, q->cell.ssb.Lmax, n_hf);
			peak_value = srslte_vec_dot_prod_conj_ccc(pilot_recv_signal, pilot_refs_signal[n_hf][i_ssb],
			SRSLTE_PBCH_DMRS_LEN);
			if (cabsf(peak_value) > cabsf(peak_max)) {
				hf_pos = n_hf;
				ssb_pos = i_ssb;
				peak_max = peak_value;
			}
		}
	}

	q->cell.ssb.half_frame = hf_pos;
	q->cell.ssb.ssb_index = ssb_pos;

	for (n_hf = 0; n_hf < 2; n_hf++) {
		memcpy(q->refs.pilots[0][n_hf][1], pilot_refs_signal[n_hf][ssb_pos] + 0, 60 * sizeof(cf_t));
		memcpy(q->refs.pilots[0][n_hf][2], pilot_refs_signal[n_hf][ssb_pos] + 60, 12 * sizeof(cf_t));
		memcpy(q->refs.pilots[0][n_hf][2] + 48, pilot_refs_signal[n_hf][ssb_pos] + 72, 12 * sizeof(cf_t));
		memcpy(q->refs.pilots[0][n_hf][3], pilot_refs_signal[n_hf][ssb_pos] + 84, 60 * sizeof(cf_t));
	}
	srslte_sequence_free(&seq);

	return cabsf(peak_max);
}

/** Initializes the PBCH transmitter and receiver.
 * At the receiver, the field nof_ports in the cell structure indicates the
 * maximum number of BS transmitter ports to look for.
 */
int srslte_pbch_init(srslte_pbch_t* q) {
	int ret = SRSLTE_ERROR_INVALID_INPUTS;

	if (q != NULL) {
		ret = SRSLTE_ERROR;

		bzero(q, sizeof(srslte_pbch_t));

		if (srslte_pbch_res_init(&q->chest.res, &q->chest.cell)) {
			goto clean;
		}
		if (srslte_modem_table_lte(&q->mod, SRSLTE_MOD_QPSK)) {
			goto clean;
		}
		if (srslte_polar_decode_init(&q->decoder)) {
			goto clean;
		}
		if (srslte_crc_init(&q->crc, SRSLTE_CRC24C, 24)) {
			goto clean;
		}

		q->nof_symbols = SRSLTE_PBCH_RE_CP_NORM;
		q->llr = srslte_vec_f_malloc(q->nof_symbols * q->mod.nbits_x_symbol);
		if (!q->llr) {
			goto clean;
		}

		ret = SRSLTE_SUCCESS;
	}
	clean: if (ret == SRSLTE_ERROR) {
		srslte_pbch_free(q);
	}
	return ret;
}

void srslte_pbch_free(srslte_pbch_t* q) {
	srslte_sequence_free(&q->seq);
	srslte_modem_table_free(&q->mod);
	srslte_polar_decode_free(&q->decoder);
	if (q->llr) {
		free(q->llr);
	}
	bzero(q, sizeof(srslte_pbch_t));
}

int srslte_pbch_scrambling_seq(srslte_sequence_t* seq, uint32_t cell_id, uint8_t ssb_index, uint8_t Lmax) {
	uint32_t v, M;

	if (Lmax == 4) {
		v = ssb_index & 0x03;
	} else {
		v = ssb_index & 0x07;
	}
	M = SRSLTE_PBCH_ENCODED_LEN;

	return srslte_sequence_LTE_pr_offset(seq, v * M, SRSLTE_PBCH_ENCODED_LEN, cell_id);
}

int srslte_pbch_payload_seq(srslte_sequence_t* seq, uint32_t cell_id, uint8_t v, uint8_t M) {
	return srslte_sequence_LTE_pr_offset(seq, v * M, SRSLTE_PBCH_PAYLOAD_LEN, cell_id);
}

int srslte_pbch_set_cell(srslte_pbch_t* q, srslte_cell_t* cell) {
	int ret = SRSLTE_ERROR_INVALID_INPUTS;

	if (q != NULL && srslte_cell_isvalid(cell)) {
		if (q->chest.cell.id != cell->id) {
			q->chest.cell = *cell;
			if (srslte_pbch_res_init(&q->chest.res, cell)) {
				return SRSLTE_ERROR;
			}
			if (srslte_chest_dl_init(&q->chest, &q->chest.res, cell)) {
				return SRSLTE_ERROR;
			}
			if (srslte_chest_cfg_init(&q->chest_cfg, 1, 1)) {
				return SRSLTE_ERROR;
			}
			if (srslte_pbch_scrambling_seq(&q->seq, cell->id, cell->ssb.ssb_index, cell->ssb.Lmax)) {
				return SRSLTE_ERROR;
			}
		}
		q->nof_symbols = SRSLTE_PBCH_RE_NUM;

		ret = SRSLTE_SUCCESS;
	}
	return ret;
}

int srslte_ssb_start_symbol(srslte_cell_t* cell) {
	srslte_cp_t cp = cell->cp;
	srslte_scs_type_t mu = cell->mu;
	uint8_t i_ssb = cell->ssb.ssb_index;
	uint8_t half_frame = cell->ssb.half_frame;
	srslte_ssb_type_t type = cell->ssb.ssb_type;
	uint32_t case_AC[2] = { 2, 8 };
	uint32_t case_BD[4] = { 4, 8, 16, 20 };
	uint32_t case_E[8] = { 8, 12, 16, 20, 32, 36, 40, 44 };
	uint32_t n, n_temp;
	uint32_t symbol;

	switch (mu) {
	case SRSLTE_MU_0: // case A
		n = i_ssb >> 1;
		symbol = case_AC[i_ssb % 2] + 14 * n;
		break;
	case SRSLTE_MU_1:
		if (type == 1) { // case B
			n = i_ssb >> 2;
			symbol = case_BD[i_ssb % 4] + 28 * n;
		}
		if (type == 2) { // case C
			n = i_ssb >> 1;
			symbol = case_AC[i_ssb % 2] + 14 * n;
		}
		break;
	case SRSLTE_MU_3: // case D
		n_temp = i_ssb >> 2;
		n = n_temp + (n_temp >> 2);
		symbol = case_BD[i_ssb % 4] + 28 * n;
		break;
	case SRSLTE_MU_4:  // case E
		n_temp = i_ssb >> 3;
		n = n_temp + (n_temp >> 2);
		symbol = case_E[i_ssb % 8] + 56 * n;
		break;
	default:
		ERROR("Error in PBCH start symbol: Invalid numerology index %d for the synchronization block\n", mu);
	}

	if (half_frame) {
		symbol += 5 * SRSLTE_CP_NSYMB(cp) * SRSLTE_N_SUBFRAME_SLOT(mu);
	}

	return symbol;
}

/**
 * Unpacks MIB from PBCH message.
 *
 * @param[in] msg PBCH in an unpacked bit array of size 24
 * @param[out] sfn System frame number
 * @param[out] cell MIB information about PHICH and system bandwidth will be saved here
 */
void srslte_pbch_mib_unpack(uint8_t* msg, srslte_cell_t* cell, uint32_t* sfn) {
	uint32_t n_f;
	uint32_t value;

	value = srslte_bit_pack(&msg, 1);
	value = srslte_bit_pack(&msg, 6);
	n_f = value << 0x04;
	value = srslte_bit_pack(&msg, 1);
	value = srslte_bit_pack(&msg, 4);
	value = srslte_bit_pack(&msg, 1);
	value = srslte_bit_pack(&msg, 8);
	value = srslte_bit_pack(&msg, 1);
	value = srslte_bit_pack(&msg, 1);
	value = srslte_bit_pack(&msg, 1);

	value = srslte_bit_pack(&msg, 4);
	n_f |= value;
	value = srslte_bit_pack(&msg, 1);
	cell->ssb.half_frame = value;
	if (cell->ssb.Lmax == 64) {
		value = srslte_bit_pack(&msg, 3);
		cell->ssb.ssb_index |= value << 0x03;
	} else {
		value = srslte_bit_pack(&msg, 1);
		msg += 2;
	}
	cell->ssb.symbol_offset = srslte_ssb_start_symbol(cell);

	if (sfn) {
		*sfn = n_f;
	}
}

void srslte_scrambling_pbch_mask(srslte_sequence_t* s, uint8_t* data, uint32_t mask, int offset, int len) {
	for (int i = 0, j = 0; i < len; i++) {
		if (!SRSLTE_BIT_GET(mask, i)) {
			data[i] ^= s->c[offset + j];
			j++;
		}
	}
}

/* Checks CRC after applying the mask for the given number of ports.
 *
 * The bits buffer size must be at least 40 bytes.
 *
 * Returns 0 if the data is correct, -1 otherwise
 */
uint32_t srslte_pbch_crc_check(srslte_pbch_t* q, uint8_t* bits, uint32_t nof_ports) {
	int ret = srslte_crc_checksum(&q->crc, bits, SRSLTE_PBCH_PAYLOADCRC_LEN);
	if (ret == 0) {

	} else {
		return ret;
	}
	return ret;
}

void srslte_pbch_deinterleaver(srslte_pbch_t* q) {
	uint8_t pattern[32] = { 28, 0, 31, 30, 7, 29, 25, 27, 5, 8, 24, 9, 10, 11, 12, 13, 1, 4, 3, 14, 15, 16, 17, 2, 26, 18,
			19, 20, 21, 22, 6, 23 };

	for (int i = 0; i < SRSLTE_PBCH_PAYLOAD_LEN; i++) {
		q->data[pattern[i]] = q->data_dec[i];
	}
}

int pbch_decode_data(srslte_pbch_t* q, uint32_t nof_bits) {
	srslte_cell_t* cell = &q->chest.cell;
	uint32_t ret;
	uint32_t Lmax = cell->ssb.Lmax;
	uint32_t ssb_index = cell->ssb.ssb_index;
	uint32_t v = Lmax == 4 ? ssb_index & 0x03 : ssb_index & 0x07;
	uint32_t M = nof_bits;
	uint32_t offset = v * M;
	uint32_t mask = Lmax == 64 ? 0x100006D : 0x1000041;

	if (srslte_pbch_scrambling_seq(&q->seq, cell->id, cell->ssb.ssb_index, cell->ssb.Lmax)) {
		return SRSLTE_ERROR;
	}
	srslte_scrambling_f_offset(&q->seq, q->llr, 0, nof_bits);

	ret = srslte_polar_set_decode(&q->decoder, 1, SRSLTE_PBCH_PAYLOADCRC_LEN, nof_bits, 0, 0, 0, 0, 0);
	ret = srslte_polar_decode(&q->decoder, q->llr, nof_bits, q->data_dec, SRSLTE_PBCH_PAYLOAD_LEN, NULL);
	if (ret == SRSLTE_SUCCESS) {
		v = q->data_dec[24] | (q->data_dec[6] << 0x01);
		M = Lmax == 64 ? SRSLTE_PBCH_PAYLOAD_LEN - 6 : SRSLTE_PBCH_PAYLOAD_LEN - 3;
		offset = v * M;
		srslte_pbch_payload_seq(&q->seq_payload, cell->id, v, M);
		srslte_scrambling_pbch_mask(&q->seq_payload, q->data_dec, mask, 0, SRSLTE_PBCH_PAYLOAD_LEN);
		srslte_pbch_deinterleaver(q);
	}
	return ret;
}

/* Decodes the PBCH channel
 *
 * The PBCH spans in 40 ms. This function is called every 10 ms. It tries to decode the MIB
 * given the symbols of a subframe (1 ms). Successive calls will use more subframes
 * to help the decoding process.
 *
 * Returns 1 if successfully decoded MIB, 0 if not and -1 on error
 */
int srslte_pbch_decode(srslte_pbch_t* q, srslte_dl_sf_cfg_t* sf, cf_t* sf_symbols[SRSLTE_MAX_PORTS],
		uint8_t bch_payload[SRSLTE_PBCH_PAYLOAD_LEN]) {
	int ret = SRSLTE_ERROR_INVALID_INPUTS;

	if (q != NULL && sf_symbols != NULL) {
		srslte_pbch_refs_corr(&q->chest, sf_symbols);
		srslte_chest_pilot_get(&q->chest, q->chest.cell.ssb.half_frame, 0);
		srslte_chest_dl_estimate(&q->chest, &q->chest_cfg, sf, sf_symbols, &q->chest.meas);
		srslte_predecoder_demap(&q->chest, &q->chest_cfg, sf_symbols);

		for (int i = 0; i < q->nof_symbols; i++) {
			printf("d: %g %g\n", __real__ q->chest.merge.d[0][i], __imag__ q->chest.merge.d[0][i]);
		}
		srslte_demod_soft_demodulate(SRSLTE_MOD_QPSK, q->chest.merge.d[0], q->llr, q->nof_symbols);
		int nof_bits = 2 * q->nof_symbols;
		for (int i = 0; i < nof_bits; i++) {
			printf("llr: %g\n", q->llr[i]);
		}

		ret = pbch_decode_data(q, nof_bits);
		if (ret == SRSLTE_SUCCESS) {
			if (bch_payload) {
				memcpy(bch_payload, q->data, sizeof(uint8_t) * SRSLTE_PBCH_PAYLOAD_LEN);
			}
		}
	}

	return ret;
}
