/*
 * pbch.c
 *
 *  Created on: Jun 25, 2020
 *      Author: administrator
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "srslte/phy/dft/dft.h"
#include "srslte/phy/sync/sss.h"
#include "srslte/phy/utils/convolution.h"
#include "srslte/phy/utils/debug.h"
#include "srslte/phy/utils/vector.h"

int srslte_sss_init(srslte_sss_t* q, uint32_t fft_size) {

	if (q != NULL && fft_size <= 2048) {
		uint32_t N_id_2;

		bzero(q, sizeof(srslte_sss_t));

		if (srslte_dft_plan(&q->dftp_input, fft_size, SRSLTE_DFT_FORWARD, SRSLTE_DFT_COMPLEX)) {
			srslte_sss_free(q);
			return SRSLTE_ERROR;
		}
		srslte_dft_plan_set_mirror(&q->dftp_input, true);
		srslte_dft_plan_set_dc(&q->dftp_input, true);

		q->fft_size = fft_size;
		q->max_fft_size = fft_size;

		for (N_id_2 = 0; N_id_2 < 3; N_id_2++) {
			generate_sss_all_tables(&q->tables[N_id_2], N_id_2);
		}
		q->N_id_2 = 0;
		return SRSLTE_SUCCESS;
	}
	return SRSLTE_ERROR_INVALID_INPUTS;
}

int srslte_sss_resize(srslte_sss_t* q, uint32_t fft_size) {
	if (q != NULL && fft_size <= 2048) {
		if (fft_size > q->max_fft_size) {
			ERROR("Error in sss_synch_resize(): fft_size must be lower than initialized\n");
			return SRSLTE_ERROR;
		}
		if (srslte_dft_replan(&q->dftp_input, fft_size)) {
			srslte_sss_free(q);
			return SRSLTE_ERROR;
		}
		q->fft_size = fft_size;
		return SRSLTE_SUCCESS;
	}
	return SRSLTE_ERROR_INVALID_INPUTS;
}

void srslte_sss_free(srslte_sss_t* q) {
	srslte_dft_plan_free(&q->dftp_input);
	bzero(q, sizeof(srslte_sss_t));
}

/** Sets the N_id_2 to search for */
int srslte_sss_set_N_id_2(srslte_sss_t* q, uint32_t N_id_2) {
	if (!srslte_N_id_2_isvalid(N_id_2)) {
		ERROR("Invalid N_id_2 %d\n", N_id_2);
		return SRSLTE_ERROR;
	} else {
		q->N_id_2 = N_id_2;
		return SRSLTE_SUCCESS;
	}
}

void srslte_sss_generate(float* signal, uint32_t N_id_2, uint32_t N_id_1) {
	int16_t m0, m1;

	int16_t x0[SRSLTE_SSS_LEN];
	int16_t x1[SRSLTE_SSS_LEN];

	int16_t x0_initial[SRSLTE_SSS_INITIAL] = { 1, 0, 0, 0, 0, 0, 0 };
	int16_t x1_initial[SRSLTE_SSS_INITIAL] = { 1, 0, 0, 0, 0, 0, 0 };

	for (int i = 0; i < SRSLTE_SSS_INITIAL; i++) {
		x0[i] = x0_initial[i];
		x1[i] = x1_initial[i];
	}

	for (int i = 0; i < (SRSLTE_SSS_LEN - SRSLTE_SSS_INITIAL); i++) {
		x0[i + 7] = (x0[i + 4] + x0[i]) & 0x01;
		x1[i + 7] = (x1[i + 1] + x1[i]) & 0x01;
	}

	m0 = 15 * (N_id_1 / 112) + (5 * N_id_2);
	m1 = N_id_1 % 112;

	for (int n = 0; n < SRSLTE_SSS_LEN; n++) {
		signal[n] = (1 - 2 * x0[(n + m0) % SRSLTE_SSS_LEN]) * (1 - 2 * x1[(n + m1) % SRSLTE_SSS_LEN]);
	}
}

void generate_sss_all_tables(srslte_sss_tables_t* tables, uint32_t N_id_2) {
	uint32_t i;

	for (i = 0; i < SRSLTE_SSS_NUM; i++) {
		srslte_sss_generate(tables->s[i], N_id_2, i);
	}
}

/** 36.211 10.3 section 6.11.2.2
 */
void srslte_sss_put_slot(float* sss, cf_t* slot, uint32_t offset, srslte_cp_t cp) {
	uint32_t i, k = offset;

	if (k > 8) {
		memset(&slot[k - 8], 0, 8 * sizeof(cf_t));
		for (i = 0; i < SRSLTE_SSS_LEN; i++) {
			__real__ slot[k + i] = sss[i];
			__imag__ slot[k + i] = 0;
		}
		memset(&slot[k + SRSLTE_SSS_LEN], 0, 9 * sizeof(cf_t));
	}
}

void srslte_sss_get_slot(cf_t* sss, cf_t* slot, uint32_t ssb_offset, srslte_cp_t cp) {
	uint32_t k;

	k = ssb_offset + SRSLTE_SSS_START;

	memcpy(sss, &slot[k], SRSLTE_SSS_LEN * sizeof(cf_t));
}

/** Sets the SSS correlation peak detection threshold */
void srslte_sss_set_threshold(srslte_sss_t* q, float threshold) {
	q->corr_peak_threshold = threshold;
}

/** Returns the subframe index based on the m0 and m1 values */
uint32_t srslte_sss_subframe(uint32_t m0, uint32_t m1) {
	return 0;
}

/** Returns the N_id_1 value based on the m0 and m1 values */
int srslte_sss_N_id_1(srslte_sss_t* q, uint32_t id, float corr) {
	int N_id_1 = SRSLTE_ERROR;

	// Check threshold, consider not found (error) if the correlation is not above the threshold
	if (corr > q->corr_peak_threshold) {
		N_id_1 = id;
	}

	return N_id_1;
}
