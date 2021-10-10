/*
 * srs.c
 *
 *  Created on: Jun 25, 2020
 *      Author: administrator
 */

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "srslte/phy/ch_estimation/refsignal_ul.h"
#include "srslte/phy/common/phy_common.h"
#include "srslte/phy/common/sequence.h"
#include "srslte/phy/dft/dft_precoding.h"
#include "srslte/phy/phch/srs.h"
#include "srslte/phy/utils/debug.h"
#include "srslte/phy/utils/vector.h"
#include "srs_tables.h"

const uint16_t SRS_antenna_port[] = { 1000, 1001, 1002, 1003 };
const uint16_t srs_period[] = { 1, 2, 4, 5, 8, 10, 16, 20, 32, 40, 64, 80, 160, 320, 640, 1280, 2560 };

int srslte_refsignal_srs_init(srslte_refsignal_srs_t* q, uint32_t max_prb) {
	return SRSLTE_SUCCESS;
}

int srslte_refsignal_srs_set_cell(srslte_refsignal_srs_t* q, srslte_cell_t* cell) {
	return SRSLTE_SUCCESS;
}

void srslte_refsignal_srs_free(srslte_refsignal_srs_t* q) {

}

static uint32_t srs_cs_max(srslte_refsignal_srs_cfg_t* cfg) {
	uint8_t n_SRS_cs_max;
	uint8_t n_SRS_cs = cfg->cyclicShift;
	uint8_t K_TC = cfg->transmissionComb;

	/* see 38211 6.4.1.4.2 Sequence generation */
	if (K_TC == 4) {
		n_SRS_cs_max = 12;
		// delta = 2;     /* delta = log2(K_TC) */
	} else if (K_TC == 2) {
		n_SRS_cs_max = 8;
		// delta = 1;     /* delta = log2(K_TC) */
	} else {
		ERROR("generate_srs: SRS unknown value for K_TC %d !\n", K_TC);
		return SRSLTE_ERROR;
	}
	if (n_SRS_cs >= n_SRS_cs_max) {
		ERROR("generate_srs: inconsistent parameter n_SRS_cs %d >=  n_SRS_cs_max %d !\n", n_SRS_cs, n_SRS_cs_max);
		return SRSLTE_ERROR;
	}

	return n_SRS_cs_max;
}

static float srs_alpha(srslte_refsignal_srs_cfg_t* cfg, uint8_t p_index) {
	float alpha_i;
	uint16_t n_SRS_cs_i;
	uint8_t n_SRS_cs_max;
	uint8_t n_SRS_cs = cfg->cyclicShift;
	uint8_t K_TC = cfg->transmissionComb;
	uint8_t N_ap = (uint8_t) cfg->nrof_SrsPorts;

	n_SRS_cs_max = srs_cs_max(cfg);
	n_SRS_cs_i = (n_SRS_cs + (n_SRS_cs_max * (SRS_antenna_port[p_index] - 1000) / N_ap)) % n_SRS_cs_max;
	alpha_i = 2 * M_PI * ((float) n_SRS_cs_i / (float) n_SRS_cs_max);

	return alpha_i;
}

static uint32_t srs_hopping(srslte_refsignal_srs_t* q, srslte_refsignal_srs_cfg_t* cfg, uint32_t n_sf, uint32_t l,
		uint32_t M_sc_b_SRS, uint32_t* u, uint32_t* v) {
	uint32_t f_gh, v_nu;
	uint8_t groupOrSequenceHopping = cfg->groupOrSequenceHopping;
	uint8_t n_ID_SRS = cfg->sequenceId;

	if (groupOrSequenceHopping != neitherHopping) {
		ERROR("generate_srs: sequence hopping is not yet supported!\n");
		return SRSLTE_ERROR;
	}

	switch (groupOrSequenceHopping) {
	case neitherHopping: {
		f_gh = 0;
		v_nu = 0;
		break;
	}
	case groupHopping: {
		f_gh = q->f_gh[n_sf][l] % 30;
		v_nu = 0;
		break;
	}
	case sequenceHopping: {
		f_gh = 0;
		if (M_sc_b_SRS > 6 * SRSLTE_NRE) {
			v_nu = q->v[n_sf][l];
		} else {
			v_nu = 0;
		}
		break;
	}
	default: {
		ERROR("generate_srs: unknown hopping setting %d !\n", groupOrSequenceHopping);
		return SRSLTE_ERROR;
	}
	}

	/* u is the sequence-group number defined in section 6.4.1.4.1 */
	*u = (f_gh + n_ID_SRS) % SRSLTE_U_GROUP_NUMBER; /* 30 */
	*v = v_nu;

	return SRSLTE_SUCCESS;
}

static uint32_t srs_Fb(srslte_refsignal_srs_cfg_t* cfg, srslte_ul_sf_cfg_t *sf, uint16_t* n_b, uint8_t lprime,
		resourceType_t res_type) {
	uint16_t N_b, F_b, n_SRS;
	/* get parameters from SRS resource configuration */
	uint8_t B_SRS = cfg->freqHopping_b_SRS;
	uint8_t C_SRS = cfg->freqHopping_c_SRS;
	uint8_t b_hop = cfg->freqHopping_b_hop;
	uint8_t n_RRC = cfg->freqDomainPosition;
	uint8_t R = cfg->resourceMapping_repetitionFactor;
	uint8_t N_symb_SRS = cfg->resourceMapping_nrofSymbols; /* consecutive OFDM symbols */
	uint16_t T_SRS = srs_period[cfg->SRS_Periodicity];
	uint16_t T_offset = cfg->SRS_Offset; /* FFS_TODO_NR to check interface with RRC */
	uint16_t m_SRS_b = srs_bandwidth_config[C_SRS][B_SRS][0]; /* m_SRS_b is given by TS 38211 clause 6.4.1.4.3 */

	if (R == 0) {
		ERROR("generate_srs: this parameter repetition factor %d is not consistent !\n", R);
		return SRSLTE_ERROR;
	} else if (R > N_symb_SRS) {
		ERROR("generate_srs: R %d can not be greater than N_symb_SRS %d !\n", R, N_symb_SRS);
		return SRSLTE_ERROR;
	}

	if (T_SRS == 0) {
		ERROR("generate_srs: inconsistent parameter T_SRS %d can not be equal to zero !\n", T_SRS);
		return SRSLTE_ERROR;
	} else {
		int index = 0;
		while (srs_periodicity[index] != T_SRS) {
			index++;
			if (index == SRS_PERIODICITY) {
				ERROR("generate_srs: inconsistent parameter T_SRS %d not specified !\n", T_SRS);
				return SRSLTE_ERROR;
			}
		}
	}

	for (int b = 0; b <= B_SRS; b++) {
		N_b = srs_bandwidth_config[C_SRS][b][1];
		if (b_hop >= B_SRS) {
			n_b[b] = (4 * n_RRC / m_SRS_b) % N_b; /* frequency hopping is disabled and the frequency position index n_b remains constant */
		} else {
			if (b == b_hop) {
				N_b = 1;
			}
			/* periodicity and offset */
			if (res_type == aperiodic) {
				n_SRS = lprime / R;
			} else {
				int8_t N_slot_frame = SRSLTE_N_FRAME_SLOT(sf->mu);
				n_SRS = ((N_slot_frame * sf->n_f + sf->n_sf - T_offset) / T_SRS) * (N_symb_SRS / R) + (lprime / R);
			}

			uint16_t product_N_b = 1; /* for b_hop to b-1 */
			for (unsigned int b_prime = b_hop; b_prime < B_SRS; b_prime++) { /* product for b_hop to b-1 */
				if (b_prime != b_hop) {
					product_N_b *= srs_bandwidth_config[C_SRS][b_prime][1];
				}
			}

			if (N_b & 1) { /* Nb odd */
				F_b = (N_b / 2) * (n_SRS / product_N_b);
			} else { /* Nb even */
				uint16_t product_N_b_B_SRS = product_N_b;
				product_N_b_B_SRS *= srs_bandwidth_config[C_SRS][B_SRS][1]; /* product for b_hop to b */
				F_b = (N_b / 2) * ((n_SRS % product_N_b_B_SRS) / product_N_b) + ((n_SRS % product_N_b_B_SRS) / 2 * product_N_b);
			}

			if (b <= b_hop) {
				n_b[b] = (4 * n_RRC / m_SRS_b) % N_b;
			} else {
				n_b[b] = (F_b + (4 * n_RRC / m_SRS_b)) % N_b;
			}
		}
	}
	return F_b;
}

static uint32_t srs_k0_ue(srslte_refsignal_srs_cfg_t* cfg, uint16_t* n_b, uint8_t p_index) {
	uint8_t n_SRS_cs_max;
	uint8_t B_SRS = cfg->freqHopping_b_SRS;
	uint8_t C_SRS = cfg->freqHopping_c_SRS;
	uint16_t n_shift = cfg->freqDomainShift;
	uint8_t n_SRS_cs = cfg->cyclicShift;
	uint8_t K_TC = cfg->transmissionComb;
	uint8_t K_TC_overbar = cfg->combOffset; /* FFS_TODO_NR is this parameter for K_TC_overbar ?? */
	uint8_t N_ap = (uint8_t) cfg->nrof_SrsPorts; /* antenna port for transmission */
	uint16_t m_SRS_b = srs_bandwidth_config[C_SRS][B_SRS][0]; /* m_SRS_b is given by TS 38211 clause 6.4.1.4.3 */
	uint16_t M_sc_b_SRS = m_SRS_b * SRSLTE_NRE / K_TC; /* length of the sounding reference signal sequence */
	uint16_t K_TC_p;
	uint16_t k_0_overbar_p;
	uint16_t k_0_p;

	n_SRS_cs_max = srs_cs_max(cfg);
	/* TS 38.211 6.4.1.4.3  Mapping to physical resources */

	if ((n_SRS_cs >= n_SRS_cs_max / 2) && (n_SRS_cs < n_SRS_cs_max) && (N_ap == 4)
			&& ((SRS_antenna_port[p_index] == 1001) || (SRS_antenna_port[p_index] == 1003))) {
		K_TC_p = (K_TC_overbar + K_TC / 2) % K_TC;
	} else {
		K_TC_p = K_TC_overbar;
	}

	k_0_overbar_p = n_shift * SRSLTE_NRE + K_TC_p;

	k_0_p = k_0_overbar_p; /* it is the frequency-domain starting position */

	for (int b = 0; b <= B_SRS; b++) {
		k_0_p += K_TC * M_sc_b_SRS * n_b[b];
	}

	return k_0_p;
}

uint32_t srslte_refsignal_srs_M_sc(srslte_refsignal_srs_cfg_t* cfg) {
	uint8_t B_SRS = cfg->freqHopping_b_SRS;
	uint8_t C_SRS = cfg->freqHopping_c_SRS;
	uint8_t K_TC = cfg->transmissionComb;

	uint16_t m_SRS_b = srs_bandwidth_config[C_SRS][B_SRS][0]; /* m_SRS_b is given by TS 38211 clause 6.4.1.4.3 */
	uint16_t M_sc_b_SRS = m_SRS_b * SRSLTE_NRE / K_TC; /* length of the sounding reference signal sequence */

	return M_sc_b_SRS;
}

int32_t generate_srs_nr(srslte_refsignal_srs_t* srs, srslte_refsignal_srs_cfg_t* cfg, srslte_ul_sf_cfg_t *sf,
		cf_t* r_srs) {
	float alpha;
	uint32_t u, v, k_0_p;
	uint16_t n_b[B_SRS_NUMBER];

	/* TS 38.211 6.4.1.4.1 SRS resource */
	uint8_t N_ap = (uint8_t) cfg->nrof_SrsPorts; /* antenna port for transmission */
	uint8_t N_symb_SRS = cfg->resourceMapping_nrofSymbols; /* consecutive OFDM symbols */
	uint8_t l_offset = cfg->resourceMapping_startPosition;
	uint8_t l0 = SRSLTE_N_SLOT_SYMB - 1 - l_offset; /* starting position in the time domain *//* frequency domain starting position */
	resourceType_t res_type = cfg->resourceType;

	if (res_type != periodic) {
		ERROR("generate_srs: only SRS periodic is supported up to now!\n");
		return SRSLTE_ERROR;
	}

	if (N_ap != port1) {
		ERROR("generate_srs: this number of antenna ports %d is not yet supported!\n", N_ap);
		return SRSLTE_ERROR;
	}

	if (N_symb_SRS != 1) {
		ERROR("generate_srs: this number of srs symbol  %d is not yet supported!\n", N_symb_SRS);
		return SRSLTE_ERROR;
	}

	uint16_t M_sc_b_SRS = srslte_refsignal_srs_M_sc(cfg);

	/* for each antenna ports for transmission */
	for (int p_index = 0; p_index < N_ap; p_index++) {
		/* see TS 38.211 6.4.1.4.2 Sequence generation */
		alpha = srs_alpha(cfg, p_index);

		/* for each SRS symbol which should be managed by SRS configuration */
		/* from TS 38.214 6.2.1.1 UE SRS frequency hopping procedure */
		/* A UE may be configured to transmit an SRS resource on  adjacent symbols within the last six symbols of a slot, */
		/* where all antenna ports of the SRS resource are mapped to each symbol of the resource */
		uint8_t l = p_index;
		if (l >= N_symb_SRS) {
			ERROR("generate_srs: number of antenna ports %d and number of srs symbols %d are different !\n", N_ap,
					N_symb_SRS);
		}

		srs_hopping(srs, cfg, sf->n_sf, l0 + l, M_sc_b_SRS, &u, &v);
		srs_Fb(cfg, sf, n_b, l, res_type);
		k_0_p = srs_k0_ue(cfg, n_b, p_index);

		srslte_r_uv_seq(&r_srs[p_index * M_sc_b_SRS], u, v, alpha, M_sc_b_SRS);
	}

	return k_0_p;
}

