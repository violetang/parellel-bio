/**
 * @file accel_em.c
 * @author Arun Sethuraman
 * @author Karin S. Dorman
 * @date Sun Jan  4 22:33:34 CST 2015
 * 
 * This file contains functions that will perform accelerated EM.
 */

#include "multiclust.h"

double step_size(options *opt, data *dat, model *mod);

/**
 * Accelerated update. TODO: add comments
 *
 * @param dat data object
 * @param mod mod object
 * @param opt options object
 * @return log likelihood
 */
double accelerated_update(options *opt, data *dat, model *mod)
{
	int i, l, k, m, m_start;
	int n_adjust = 0;
	double ppll, pll, ll=0, s;

	/* three steps of EM: init_, iter1_, iter2_ */
	if (opt->fresh_em)
		ppll = em_3_steps(mod, dat, opt, mod->logL);
	else {
		ppll = em_1_step(mod, dat, opt, mod->logL);
		mod->accel_abort = 0;
	}
	pll = log_likelihood(opt, dat, mod, 3);

	/* calculate step size */
	s = step_size(opt, dat, mod);
//fprintf(stdout, "Step size is %f\n", s);

	/* fall back on EM step */
	if (isnan(s))
		goto EM_EXIT;

/*
print_param(opt, dat, mod, 1);
print_param(opt, dat, mod, 2);
print_param(opt, dat, mod, 3);
*/

	do {

		/* update allele frequencies */
		for (l = 0; l < dat->L; l++)
			for (k = 0; k < mod->K; k++) {
				m_start = dat->L_alleles
					&& dat->L_alleles[l][0] == MISSING;
				for (m=m_start; m < dat->uniquealleles[l]; m++)
					if (opt->accel_scheme == QN1)
						mod->pKLM[k][l][m] = 
							(1-s) * mod->iter1_pKLM[k][l][m]
							+ s * mod->iter2_pKLM[k][l][m];
					else
						mod->pKLM[k][l][m] = 
							mod->init_pKLM[k][l][m]
							- 2*s*(mod->iter1_pKLM[k][l][m] - mod->init_pKLM[k][l][m])
							+ s*s*(mod->iter2_pKLM[k][l][m] - 2*mod->iter1_pKLM[k][l][m] + mod->init_pKLM[k][l][m]);
				simplex_project_pklm(mod, dat, opt, k, l);
			}
		
		/* update mixing proportions */
		if (opt->admixture && !opt->eta_constrained) {
			for (i = 0; i < dat->I; i++) {
				for (k = 0; k < mod->K; k++)
					if (opt->accel_scheme == QN1)
						mod->etaik[i][k] =
							(1-s) * mod->iter1_etaik[i][k]
							+ s * mod->iter2_etaik[i][k];
					else
						mod->etaik[i][k] =
							mod->init_etaik[i][k]
							- 2 * s * (mod->iter1_etaik[i][k] - mod->init_etaik[i][k])
							+ s * s * (mod->iter2_etaik[i][k] - 2*mod->iter1_etaik[i][k] + mod->init_etaik[i][k]);
				simplex_project_eta(mod, opt, i);
			}
		} else {
			for (k = 0; k < mod->K; k++)
				if (opt->accel_scheme == QN1)
					mod->etak[k] =
						(1-s) * mod->iter1_etak[k]
						+ s * mod->iter2_etak[k];
				else
					mod->etak[k] =
						mod->init_etak[k]
						- 2 * s * (mod->iter1_etak[k] - mod->init_etak[k])
						+ s * s * (mod->iter2_etak[k] - 2*mod->iter1_etak[k] + mod->init_etak[k]);
			simplex_project_eta(mod, opt, -1);
		}

//print_param(opt, dat, mod, 0);
		ll = log_likelihood(opt, dat, mod, 0);
		if (opt->verbosity)
			fprintf(stderr, "after attempt %d of accel ll is %f, EM ll is %f, step is %f\n", n_adjust, ll, pll, s);

		if (opt->adjust_step && ll < pll)
			s = (s - 1)/2;

	} while (n_adjust++ < opt->adjust_step && ll < pll && s < -1);

	/* accelerated step improved log likelihood */
	if (ll > pll) {
		mod->logL = ppll;
		if (opt->verbosity)
			fprintf(stderr, "accelerated_update (%s): %f (step size: %f; EM ll: %f)\n", 
				accel_method_abbreviations[opt->accel_scheme], ll, s, pll);
		return ll;
	}

EM_EXIT:
	if (opt->verbosity)
		fprintf(stderr, "accelerated_update (EM): %f (step size: %f; %s ll: %f)\n",
			pll, s, accel_method_abbreviations[opt->accel_scheme], ll);
	mod->accel_abort = 1;
	COPY_3JAGGED_ARRAY(mod->pKLM, mod->iter2_pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->etaik, mod->iter2_etaik, mod->K);
	else
		COPY_1ARRAY(mod->etak, mod->iter2_etak, mod->K);

	ll = pll;
	mod->logL = ppll;
	mod->n_iter--;

	return ll;
} /* End of accelerated_update */

double step_size(options *opt, data *dat, model *mod)
{
	int i, l, k, m;
	int m_start;
	double utu = 0;
	double utvu = 0;
	double vutvu = 0;
	double s;

	if (opt->admixture && !opt->eta_constrained) {
		for (i = 0; i < dat->I; i++) {
			mod->current_i = i;
			initialize_etaiks(mod);
			for (k = 0; k < mod->K; k++) {	
				utu += mod->U[k] * mod->U[k];
				utvu += mod->U[k] * (mod->V[k] - mod->U[k]);
				vutvu += (mod->V[k] - mod->U[k]) * (mod->V[k] - mod->U[k]);
			}
//fprintf(stderr, "step_size: %f, %f, %f\n", utu, utvu, vutvu);
		}
	} else {
		initialize_etaks(mod);
		for (k = 0; k < mod->K; k++) {	
			utu += mod->U[k] * mod->U[k];
			utvu += mod->U[k] * (mod->V[k] - mod->U[k]);
			vutvu += (mod->V[k] - mod->U[k]) * (mod->V[k] - mod->U[k]);
		}
	}

	for (k = 0; k < mod->K; k++) {
		for (l = 0; l < dat->L; l++) {
			mod->current_k = k;
			mod->current_l = l;
			initialize_pklas(mod, dat);
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				utu += mod->U[m] * mod->U[m];
				utvu += mod->U[m] * (mod->V[m] - mod->U[m]);
				vutvu += (mod->V[m] - mod->U[m]) * (mod->V[m] - mod->U[m]);
			}
//fprintf(stderr, "step_size: %f, %f, %f\n", utu, utvu, vutvu);
		}
	}
//fprintf(stderr, "step_size: %f, %f, %f\n", utu, utvu, vutvu);
	if (opt->accel_scheme == SQS1)
		s = utu / utvu;
	else if (opt->accel_scheme == SQS2)
		s = utvu / vutvu;
	else if (opt->accel_scheme == SQS3) {
		if (sqrt(utu) < 1e-8)
			return NAN;
		s = -sqrt( utu / vutvu );
	} else if (opt->accel_scheme == QN1)
		s = - utu / utvu;
	else s = -1;

	//if (s > -1)
	//	s = -1;

	return s;
} /* step_size */
