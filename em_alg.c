/**
 * @file em_alg.c
 * @author Wei-Chen Chen
 * @author Arun Sethuraman
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Wed Dec  5 19:05:51 CST 2012
 *
 * This file contains the functions for running iterations of the EM algorithm
 * under the (ad)mixture model.  Public functions are em(), which runs the EM to
 * convergence, and em_step(), which runs a single step.
 *
 * TODO
 * - There is code to deal with numerical errors in E step (admixture), but it
 *   is commented out for time efficiency.  We may want to check for numerical
 *   problems somehow and trigger this code when appropriate.
 */

#include "multiclust.h"

double e_step_admixture_new(options *opt, data *dat, model *mod);
double e_step_admixture_orig(options *opt, data *dat, model *mod);
void m_step_admixture_new(options *opt, data *dat, model *mod);
void m_step_admixture_orig(options *opt, data *dat, model *mod);
double e_step_mixture(data *dat, model *mod);
void m_step_mixture(options *opt, data *dat, model *mod);
double scale_log_sum(double *v, int n, double max_v);

#define e_step_admixture(A, B, C) e_step_admixture_orig(A, B, C)
#define m_step_admixture(A, B, C) m_step_admixture_orig(A, B, C)

/**
 * EM algorithm.  This function iterates EM algorithm until convergence of the
 * log likelihood or until _options::max_iter.
 *
 * @param opt options object
 * @param dat data object
 * @param mod mod object
 * @return void
 */
void em(options *opt, data *dat, model *mod)
{
	int stop = 0;
	double loglik;		/* log likelihood at next step */

	mod->logL = -Inf;
	mod->n_iter = 0;

	do {
		/* one EM iteration */
		loglik = em_step(opt, dat, mod);
		mod->n_iter++;

		if (loglik < mod->logL) {	/* DEBUGGING */
			message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
				CUSTOM_ERROR, "log likelihood decrease");
			exit(0);
		}

		/* check for convergence */
		stop = stop_iterations(opt, mod, loglik);

		if (opt->verbosity)
			fprintf(stderr, "%4d (EM): %.2f (delta): %.5g\n", mod->n_iter, loglik, loglik - mod->logL);

		/* reset */
		mod->logL = loglik;
	} while (!stop);
} /* End of em(). */


/**
 * [TODO] Comment this function.
 *
 * This is all code for block relaxation using Zhou's algorithm.
 */
void accelerated_em(options *opt, data *dat, model *mod)
{
	double loglik;
	int stop;

	/* [KSD] Simple call to em() should work; em() checks K=1 */
	if (mod->K == 1) {
		em(opt, dat, mod);
		return;
	}

	mod->n_iter = 0;
	mod->logL = -Inf;

	/* run a few EM iterations to get close */
	if (opt->n_init_iter > 3) {
		do {
			loglik = em_step(opt, dat, mod);
			mod->n_iter++;
			stop = stop_iterations(opt, mod, loglik);
			mod->logL = loglik;
		} while (mod->n_iter < opt->n_init_iter - 3 && !stop);
	}


	/* Run EM thrice to get initial estimates of all parameters */
	if (opt->block_relax)
		loglik = em_3_steps(mod, dat, opt, mod->logL);

	do {
		
//		mod->current_n = mod->n_iter % 3;

		if (opt->block_relax)
			loglik = block_update_params(opt, dat, mod); 
		else
			loglik = accelerated_update(opt, dat, mod); 

		if (isnan(loglik)) {
			message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, CUSTOM_ERROR, "nan\n");
			exit(0);
		}
		mod->n_iter++;

		stop = stop_iterations(opt, mod, loglik);

		if (opt->verbosity)
			fprintf(stderr, "%4d (%s): %.2f (delta): %.5g\n",
				mod->n_iter,
				accel_method_abbreviations[opt->accel_scheme],
				loglik, loglik - mod->logL);
//if (mod->n_iter > 21) exit(0);

		mod->logL = loglik;
		
	} while (!stop);
} /* accelerated_em */


/**
 * Check if iterations should stop.
 *
 * @param opt options object
 * @param mod model object
 * @param loglik most recent computed log likelihood; previous in mod->logL
 * @return non-zero to stop
 */
int stop_iterations(options *opt, model *mod, double loglik)
{
	double abs_diff = 0;	/* absolute difference */
	double rel_diff = 0;	/* relative difference */
	int stop = 1;

	/* K=1 converges right away (two iterates for maximum log likelihood) */
	if (mod->K == 1 && mod->n_iter == 2)
		return stop;

	if (opt->abs_error)
		abs_diff = fabs(loglik - mod->logL);
	if (opt->rel_error)
		rel_diff = abs_diff / fabs(mod->logL);

	/* [DEBUG]: verbose output for each iteration */
/*
	printf("\titer = %d, loglt = %f, loglt1 = %f, aerr = %e, rerr = %e\n",
		mod->n_iter, (float) mod->logL, loglik, abs_diff, rel_diff);
*/

	/* any one of these triggers can make the iterations continue */
	if (opt->abs_error && abs_diff > opt->abs_error)
		stop &= 0;
	if (opt->rel_error && rel_diff > opt->rel_error)
		stop &= 0;
	/* or stop for good */
	if (opt->max_iter && mod->n_iter > opt->max_iter)
		stop = 1;
	return stop;
} /* stop_iterations */

/**
 * One iteration of EM.  This function does one iteration of the EM algorithm.
 *
 * @param opt options object
 * @param dat data object
 * @param mod mod object
 * @return log likelihood from previous iteration
 */
double em_step(options *opt, data *dat, model *mod)
{
	double ll;
	if (opt->admixture) {
		ll = e_step_admixture(opt, dat, mod);
		m_step_admixture(opt, dat, mod);
	} else {
		ll = e_step_mixture(dat, mod);
		m_step_mixture(opt, dat, mod);
	}
	return ll;
} /* End of em_step(). */


/**
 * One iteration of EM followed by one E step.  The second E step is needed to
 * get the updated log likelihood following one iteration of EM.
 *
 * @param opt options object
 * @param dat data object
 * @param mod mod object
 * @return log likelihood from previous iteration
 */
double em_e_step(options *opt, data *dat, model *mod)
{
	double ll;

	if (opt->admixture) {
		/* [KSD] Warning!  Don't delete e_step_admixture() call! */
//fprintf(stderr, "first e_step: %f\n", e_step_admixture(opt, dat, mod));
		e_step_admixture(opt, dat, mod);
		m_step_admixture(opt, dat, mod);
//fprintf(stderr, "next step: %f\n", log_likelihood(opt, dat, mod, 0));
		ll = e_step_admixture(opt, dat, mod);
//fprintf(stderr, "second e_step: %f\n", ll);
	} else {
		e_step_mixture(dat, mod);
		m_step_mixture(opt, dat, mod);
		ll = e_step_mixture(dat, mod);
	}
	return ll;
} /* End of em_e_step(). */

double e_step_admixture_new(options *opt, data *dat, model *mod)
{
	int i, k, l, a, m;
	double tmp;
	double loglik = 0;

	if (!dat->ila && make_ila(dat))
		return 0;

	for (i = 0; i < dat->I; i++) {
		for (l = 0; l < dat->L; l++) {
			for (a = 0; a < dat->ploidy; a++) {
				m = dat->ila[i][l][a];

				/* [TODO] [KSD] cannot handle missing alleles */
				if (m == MISSING) {
					message(stderr, __FILE__, __func__,
						__LINE__, ERROR_MSG,
						CUSTOM_ERROR,
						"e_step_admixture() does not "
						"handle missing data");
					exit(0);
				}
				tmp = 0;
				for (k = 0; k < mod->K; k++) {
					mod->diklm[i][k][l][m] =
						(opt->eta_constrained ?
						mod->etak[k] : mod->etaik[i][k])
						* mod->pKLM[k][l][m];
					tmp += mod->diklm[i][k][l][m];
				}
				for (k = 0; k < mod->K; k++)
					mod->diklm[i][k][l][m] /= tmp;
				loglik += log(tmp);
			}
		}
	}

	return loglik;
} /* e_step_admixture_new */


/**
 * E step.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return previous step log likelihood
 */
double e_step_admixture_orig(options *opt, data *dat, model *mod)
{
	int debug = 0;
	int i, k, l, m, m_start;
	double ldilmk[mod->K];
/*	
	double scale, max_ldilmk;
*/
	double tmp, loglik = 0;
/*	
	double logeta[mod->K];
*/

	/* TIME COMPLEXITY (let T = A_1 + ... + A_L) */
	/* current version: I*K*T*(2*K+3) \propt I*K^2*T */
	/* faster version: I*K*L*M*(2*K + 2) */
	/* compare for full EM cycle):		current	vs. faster
		5*I*K*T + 2*I*K^2*T + 2*I*K + 2*K*T	vs. 2*I*K^2*L*M + 4*I*K*L*M + 2*I*K + I*K*M + 2*K*T
		5*I*K*T + 2*I*K^2*T			vs. 2*I*K^2*L*M + 4*I*K*L*M + I*K*M
		5*T + 2*K*T				vs. 2*K*L*M + 5*L*M*(4/5 + 1/(4*L))
		T*(5 + 2*K)				vs. L*M*[2*K + 5*(4/5 + 1/(4*L))]
		Suppose T = alpha*L*M, then the fast version is faster if
		alpha 					> 4/5 + 1/(4*L)
		which is almost certainly true except for SNP data on diploids
	 */
/*
	if (opt->eta_constrained)
		for (k = 0; k < mod->K; k++)
			logeta[k] = log(mod->etak[k]);
*/
	for (i = 0; i < dat->I; i++) {
/*
		if (!opt->eta_constrained)
			for (k = 0; k < mod->K; k++)
				logeta[k] = log(mod->etaik[i][k]);
*/
		for (l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING
				? 1 : 0;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				if(dat->ILM[i][l][m] == 0) {
					for (k = 0; k < mod->K; k++)
						mod->diklm[i][k][l][m] = 0;
					continue;
				}
/*
				max_ldilmk = -Inf;
*/
				tmp = 0;
				if (debug)
					fprintf(stderr, "%d, %d, %d:", i, l, m);
				for (k = 0; k < mod->K; k++) {
					ldilmk[k] = (opt->eta_constrained ?
						mod->etak[k] : mod->etaik[i][k])
						* mod->pKLM[k][l][m];
					if (debug)
						fprintf(stderr, " %f*%f",
							mod->pKLM[k][l][m],
							mod->etaik[i][k]);
/*
					ldilmk[k] = logeta[k]
						+ log(mod->pKLM[k][l][m]);
					if (ldilmk[k] > max_ldilmk)
						max_ldilmk = ldilmk[k];
					tmp += exp(ldilmk[k]);
*/
					tmp += ldilmk[k];
				}
				if (debug)
					fprintf(stderr, " = %f\n", tmp);
/*
				scale = scale_log_sum(ldilmk, mod->K, max_ldilmk);
*/
				for (k = 0; k < mod->K; k++)
					mod->diklm[i][k][l][m]
						= dat->ILM[i][l][m]
						* ldilmk[k]/tmp;
/*
					mod->diklm[i][k][l][m]
						= dat->ILM[i][l][m]
						* exp(ldilmk[k])/tmp;
*/
				loglik += dat->ILM[i][l][m] * log(tmp)
				;
/*
				+ scale;
*/
			}
			if ((!dat->L_alleles
				|| dat->L_alleles[l][0] == MISSING)
				&& dat->ILM[i][l][0] > 0) {
				for (m = m_start;
					m < dat->uniquealleles[l]; m++){
					for (k=0; k < mod->K; k++) {
						mod->diklm[i][k][l][m] += 
							(opt->eta_constrained
							? mod->etak[k]
							: mod->etaik[i][k])
							* mod->pKLM[k][l][m]
							* dat->ILM[i][l][0];
					}}

				/* TODO: FASTER?
				for (m = 0; m < M; m++) {
					temp = 0.0;
					for (j = 0; j < K; j+=)
						temp += etaik[i][j]
							* KLM[j][l][X[i][l][m]];
					diklm[i][k][l][X[i][l][m]] = etaik[i][k]
						* KLM[k][l][X[i][l][m]] / temp;
				}
				(see M-Step for matching changes) */
			}
		}
	}

	/* ALTERNATIVE: E+M combined uses less memory as I->infinity and less time as I,L->infinity
	 * Let T = A_1+...+A_L and U = max_l A_l
	 * MEMORY USAGE: NEW vs. OLD
	 *	2*I*K + 2*K*T + U*(1+K)	vs. I*K*[1 + L*M + T/I]
	 *	2*I*K*[1+T/I] + U*(1+K)	vs. I*K*[1 + L*M + T/I]
	 *	I*K*[1 + T/I] + U*(1+K)	vs. I*K*L*M
	 *	K*T + U*(1+K)		vs. I*K*(L*M-1)
	 * TIME COMPLEXITY: NEW vs. OLD
	 *	I*(1 + L*(2*calloc + 5*K*M) + 6*K*T + K) + 2*K*T	vs. I*K*T*(2*K+3) + I*K*(2+T) + K*T*(I+2)
	 *	I + I*L*(2*calloc) + 5*I*L*K*M + 6*I*K*T + I*K + 2*K*T	vs. 5*I*K*T + 2*I*K + 2*I*T*K^2 + 2*K*T
	 *	I + I*L*(2*calloc) + 5*I*L*K*M + I*K*T			vs. I*K + 2*I*T*K^2
	 */
	/*
	double *denom;
	double *dilka;
	double sum, tmp;
	int itmp;
	for (i = 0; i < I; i++) {
		sum = 0;
		for (l = 0; l < L; l++) {
			denom = calloc(uniquealleles[l], sizeof *denom);
			dilka = calloc(K*uniquealleles[l], sizeof *dilka);
			for (k = 0; k < K; k++)
				for (m = 0; m < M; m++) {
					tmp = etaik[i][k] * pkla[k][l][X[i][l][m]]
					dilka[K*X[i][l][m] + k] += tmp;
					denom[X[i][l][m]] += tmp;
				}
			for (m = 0; m < uniquealleles[l]; m++)
				for (k = 0; k < K; k++) {
					itmp = K*m + k;
					dilka[itmp] /= denom[m];
					new_pkla[k][l][m] += dilka[itmp];
					next_etaik[i][k] += dilka[itmp];
					sum += dilka[itmp];
				}
		}
		for (k = 0; k < K; k++)
			next_etaik[i][k] /= sum;
	}
	for (k = 0; k < K; k++)
		for (l = 0; l < L; l++) {
			sum = 0;
			for (m = 0; m < uniquealleles[l]; m++)
				sum += new_pkla[k][l][m];
			for (m = 0; m < uniquealleles[l]; m++)
				new_pkla[k][l][m] /= sum;
		}
	*/

	return loglik;
} /* End of e_step_admixture_orig(). */


void m_step_admixture_new(options *opt, data *dat, model *mod)
{
	int i, m, l, k, a, nhaplotypes;
	double tmp;

	for (l = 0; l < dat->L; l++) {
		if (dat->L_alleles && dat->L_alleles[l][0] == MISSING) {
			message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
				CUSTOM_ERROR, "m_step_admixture can not handle "
				"missing data");
			exit(0);
		}
		for (m = 0; m < dat->uniquealleles[l]; m++)
			for (k = 0; k < mod->K; k++)
				mod->pKLM[k][l][m] = 0;
		for (k = 0; k < mod->K; k++) {
			tmp = 0;
			for (i = 0; i < dat->I; i++)
				for (a = 0; a < dat->ploidy; a++) {
					m = dat->ila[i][l][a];
					mod->pKLM[k][l][m] +=
						mod->diklm[i][k][l][m];
					tmp += mod->diklm[i][k][l][m];
				}

			for (m = 0; m < dat->uniquealleles[l]; m++)
				mod->pKLM[k][l][m] /= tmp;
		}
	}

	if (!opt->eta_constrained) {
		nhaplotypes = dat->ploidy * dat->L;
		for (i = 0; i < dat->I; i++) {
			for (k = 0; k < mod->K; k++) {
				mod->etaik[i][k] = 0;
				for (l = 0; l < dat->L; l++)
					for (a = 0; a < dat->ploidy; a++) {
						m = dat->ila[i][l][a];
						mod->etaik[i][k] +=
							mod->diklm[i][k][l][m];
					}
				mod->etaik[i][k] /= nhaplotypes;
			}
		}
	} else {
		nhaplotypes = dat->ploidy * dat->L * dat->I;
		for (k = 0; k < mod->K; k++)
			mod->etak[k] = 0;
		for (i = 0; i < dat->I; i++)
			for (a = 0; a < dat->ploidy; a++) {
				m = dat->ila[i][l][a];
				mod->etak[k] += mod->diklm[i][k][l][m];
			}
		for (k = 0; k < mod->K; k++)
			mod->etak[k] /= nhaplotypes;
	}
} /* m_step_admixture_new */


/**
 * M step for admixture model, without accounting for missing data.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return void
 */
void m_step_admixture_orig(options *opt, data *dat, model *mod)
{
	int debug = 0;
	int i, k, l, m, m_start;
	double temp;

	/* TIME COMPLEXITY: I*K*(2+T) + K*T*(I+2) */
	/* FASTER VERSION: I*K*(2 + M) + K*(2*T + 2*L*M*I) */

	/* estimate new eta_k or eta_{ik} */
	if (!opt->admixture || opt->eta_constrained) {
		temp = 0.0;
		for (k = 0; k < mod->K; k++) {
			mod->etak[k] = 0; //opt->eta_lower_bound;
			for (l = 0; l < dat->L; l++) {
				m_start = dat->L_alleles
					&& dat->L_alleles[l][0] == MISSING
					? 1 : 0;
				for (m = m_start; m < dat->uniquealleles[l];
					m++)
					for (i = 0; i < dat->I; i++)
						mod->etak[k] +=
							mod->diklm[i][k][l][m];
			}
			temp += mod->etak[k];
		}
		if (debug)
			fprintf(stderr, "etak:");
		for (k = 0; k < mod->K; k++) {
			mod->etak[k] /= temp;
			if (debug)
				fprintf(stderr, " %f", mod->etak[k]);
		}
		if (debug)
			fprintf(stderr, "\n");
		simplex_project_eta(mod, opt, 0);
	} else {
		for (i = 0; i < dat->I; i++) {
			temp = 0.0;
			for (k = 0; k < mod->K; k++) {
				mod->etaik[i][k] = 0;//opt->eta_lower_bound;
				for (l = 0; l < dat->L; l++) {
					m_start = dat->L_alleles
						&& dat->L_alleles[l][0] == MISSING;
					for (m = m_start; m < dat->uniquealleles[l]; m++) {
						mod->etaik[i][k]
							+= mod->diklm[i][k][l][m];
						if (debug & 1)
							fprintf(stderr, "diklm[%d][%d][%d][%d]: %f (%f)\n",
								i, k, l, m, mod->diklm[i][k][l][m],
								mod->etaik[i][k]);
					}
				}
				temp += mod->etaik[i][k];
			}
			if (debug & 1)
				fprintf(stderr, "etaik[%d]:", i);
			for (k = 0; k < mod->K; k++) {
				mod->etaik[i][k] /= temp;
				if (debug & 1)
					fprintf(stderr, " %f", mod->etaik[i][k]);
			}
			if (debug & 1)
				fprintf(stderr, "\n");
			simplex_project_eta(mod, opt, i);
		}
	}

	// estimate new p_KLMs 
	for (k = 0; k < mod->K; k++) {
		for (l = 0; l < dat->L; l++) {
			temp = 0.0;
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				mod->pKLM[k][l][m] = 0;//opt->p_lower_bound;
				for (i = 0; i < dat->I; i++) {
					if (debug & 4)
						fprintf(stderr, "diklm[%d][%d][%d][%d]: %f\n", i, k, l, m, mod->diklm[i][k][l][m]);
					mod->pKLM[k][l][m] += mod->diklm[i][k][l][m];
				}
				temp += mod->pKLM[k][l][m];
			}
			if (debug & 2)
				fprintf(stderr, "pKLM[%d][%d]:", k, l);
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				mod->pKLM[k][l][m] /= temp; 
				if (debug & 2)
					fprintf(stderr, " %f", mod->pKLM[k][l][m]);
			}
			if (debug & 2)
				fprintf(stderr, "\n");
			simplex_project_pklm(mod, dat, opt, k, l);
		}
	}
} /* End of m_step_admixture */
	
/**
 * E step for mixture model.
 *
 * @param dat data object
 * @param mod model object
 * @return previous step log likelihood
 */
double e_step_mixture(data *dat, model *mod)
{
	int i, l, m, k, m_start;
	double log_etak[mod->K], temp_vik[mod->K];
	double temp;
	double scale, max_ll = -Inf;
	double loglik = 0;

	for (k = 0; k < mod->K; k++)
		log_etak[k] = log(mod->etak[k]);

	for (i = 0; i < dat->I; i++) {
		/* numerators */
		for (k = 0; k < mod->K; k++) {
			temp_vik[k] = log_etak[k];
			for (l = 0; l < dat->L; l++) {
				m_start = dat->L_alleles
					&& dat->L_alleles[l][0] == MISSING
					? 1 : 0;
				for (m = m_start; m < dat->uniquealleles[l]; m++) {
					if(dat->ILM[i][l][m] == 0)
						continue;
					temp_vik[k] += dat->ILM[i][l][m]
						* log(mod->pKLM[k][l][m]);
				}
			}
			if (temp_vik[k] > max_ll)
				max_ll = temp_vik[k];
		}

		/* normalize (possibly scaling) */
		scale = scale_log_sum(temp_vik, mod->K, max_ll);
		temp = 0;
		for(k = 0; k < mod->K; k++)
			temp += exp(temp_vik[k]);
		for(k = 0; k < mod->K; k++)
			mod->vik[i][k] = exp(temp_vik[k])/temp;

		/* restore scaling to log likelihood */
		loglik += log(temp) + scale;
	}

	return loglik;
} /* End of e_step_mixture(). */

/**
 * M step for admixture model, without accounting for missing data.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return void
 */
void m_step_mixture(options *opt, data *dat, model *mod)
{
	int i, l, m, k, m_start;
	double temp;

	/* estimate eta */
	temp = 0.0;
	for (k = 0; k < mod->K; k++) {
		mod->etak[k] = 0;//opt->eta_lower_bound;
		for (i = 0; i < dat->I; i++)
			mod->etak[k] += mod->vik[i][k];
		temp += mod->etak[k];
	}
	for (k = 0; k < mod->K; k++)
		mod->etak[k] /= temp;
	simplex_project_eta(mod, opt, 0);

	/* estimate p_{klm} */
	for (k = 0; k < mod->K; k++)
		for(l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING ? 1 : 0;
			temp = 0.0;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				mod->pKLM[k][l][m] = 0;//opt->p_lower_bound;
				for (i = 0; i < dat->I; i++) {
					if(dat->ILM[i][l][m])
						mod->pKLM[k][l][m] +=
							mod->vik[i][k]
							* dat->ILM[i][l][m];
				}
				temp += mod->pKLM[k][l][m];
			}
			for (m = m_start; m < dat->uniquealleles[l]; m++)
				mod->pKLM[k][l][m] /= temp;
			simplex_project_pklm(mod, dat, opt, k, l);
		}
} /* End of m_step_mixture(). */


/**
 * Scale vector elements.  In the EM algorithm there comes a time when
 * probabilities need to be computed by renormalizing values computed on an
 * arbitrary scale.  During normalization, we must return to the exp scale and
 * consequently there is the risk of numerical underflow (rounding to 0 of very
 * small numbers) or overflow (rounding to the biggest representable number for
 * very large numbers).  In the first case, if all numbers round to 0, then
 * there will be a divide by zero error.  In the second case, if the largest
 * number is truncated, then it will renormalize to 1.  To alleviate the
 * problem, we can translate the logged numbers so they fall within range when
 * exponentiated.  The resulting scaling factor on the non-logged scale cancels
 * during normalization.  This code finds an appropriate translation and
 * applies it so that subsequent normalization on the non-log scale will be
 * stable.
 *
 * @param v vector of logged values to normalize
 * @param n length of vector
 * @param max_v largest element of vector v
 * @return scaling parameter used
 */
double scale_log_sum(double *v, int n, double max_v)
{
	double max_ev = exp(max_v);
	double scale = 0;
	int k;

	if (max_ev == 0.0 || max_ev == HUGE_VAL) {
		/* find scaling factor that brings max_ev in range */
		scale = (max_ev == HUGE_VAL) ? max_v : -max_v;
		do {
			scale *= 0.5;
			max_ev = exp(scale);
		} while (max_ev == HUGE_VAL);

		/* diff on log scale that puts exp() in range */
		scale = max_v - scale;

		/* apply scale to vector */
		for (k=0; k<n; k++)
			v[k] -= scale;
	}
	return scale;
} /* scale_log_sum */

/**
 * Take three EM steps and store parameter estimates.  Accelerated steps
 * require storage of previous three iterations.  This function takes these
 * steps and stores the parameters.
 *
 * @param mod model object
 * @param dat data object 
 * @param opt options object
 * @return double
 */
double em_3_steps(model *mod, data *dat, options *opt, double in_ll)
{
	double ll, pll = in_ll;

	ll = em_step(opt, dat, mod);
	COPY_3JAGGED_ARRAY(mod->init_pKLM, mod->pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->init_etaik, mod->etaik, mod->K);
	else
		COPY_1ARRAY(mod->init_etak, mod->etak, mod->K);
	mod->n_iter++;

	if (opt->verbosity && (!opt->accel_scheme || mod->n_iter==1))
		fprintf(stderr, "%4d (EM): %.2f (delta): %.5g\n", mod->n_iter, ll, ll - pll);
	pll = ll;

	ll = em_step(opt, dat, mod);
	COPY_3JAGGED_ARRAY(mod->iter1_pKLM, mod->pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->iter1_etaik, mod->etaik, mod->K);
	else
		COPY_1ARRAY(mod->iter1_etak, mod->etak, mod->K);
	mod->n_iter++;

	if (opt->verbosity)
		fprintf(stderr, "%4d (EM): %.2f (delta): %.5g\n", mod->n_iter, ll, ll - pll);
	pll = ll;

	ll = em_step(opt, dat, mod);
	COPY_3JAGGED_ARRAY(mod->iter2_pKLM, mod->pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->iter2_etaik, mod->etaik, mod->K);
	else
		COPY_1ARRAY(mod->iter2_etak, mod->etak, mod->K);
	mod->n_iter++;

	if (opt->verbosity)
		fprintf(stderr, "%4d (EM): %.2f (delta): %.5g\n", mod->n_iter, ll, ll - pll);

	return ll;
} /* em_3_steps */

/**
 * Take one EM step and store parameter estimates.  Accelerated steps
 * require storage of previous three iterations.  This function takes these
 * steps and stores the parameters.
 *
 * @param mod model object
 * @param dat data object 
 * @param opt options object
 * @return double
 */
double em_1_step(model *mod, data *dat, options *opt, double in_ll)
{
	double ll, pll = in_ll;

	if (mod->accel_abort) {	/* last 3: init_, iter1_, iter2_ */
		COPY_3JAGGED_ARRAY(mod->init_pKLM, mod->iter1_pKLM, dat->uniquealleles);
		if (opt->admixture && !opt->eta_constrained)
			COPY_2ARRAY(mod->init_etaik, mod->iter1_etaik, mod->K);
		else
			COPY_1ARRAY(mod->init_etak, mod->iter1_etak, mod->K);
	
		COPY_3JAGGED_ARRAY(mod->iter1_pKLM, mod->iter2_pKLM, dat->uniquealleles);
		if (opt->admixture && !opt->eta_constrained)
			COPY_2ARRAY(mod->iter1_etaik, mod->iter2_etaik, mod->K);
		else
			COPY_1ARRAY(mod->iter1_etak, mod->iter2_etak, mod->K);
	} else {		/* last 3: iter1_, iter2_, <empty> */
		COPY_3JAGGED_ARRAY(mod->init_pKLM, mod->iter2_pKLM, dat->uniquealleles);
		if (opt->admixture && !opt->eta_constrained)
			COPY_2ARRAY(mod->init_etaik, mod->iter2_etaik, mod->K);
		else
			COPY_1ARRAY(mod->init_etak, mod->iter2_etak, mod->K);
	
		COPY_3JAGGED_ARRAY(mod->iter1_pKLM, mod->pKLM, dat->uniquealleles);
		if (opt->admixture && !opt->eta_constrained)
			COPY_2ARRAY(mod->iter1_etaik, mod->etaik, mod->K);
		else
			COPY_1ARRAY(mod->iter1_etak, mod->etak, mod->K);
	}

	ll = em_step(opt, dat, mod);
	COPY_3JAGGED_ARRAY(mod->iter2_pKLM, mod->pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->iter2_etaik, mod->etaik, mod->K);
	else
		COPY_1ARRAY(mod->iter2_etak, mod->etak, mod->K);
	mod->n_iter++;

	if (opt->verbosity && !opt->accel_scheme)
		fprintf(stderr, "%4d (EM): %.2f (delta): %.5g\n", mod->n_iter, ll, ll - pll);

	return ll;
} /* em_1_step */

/**
 * Progress forward two EM steps and store the results.
 *
 * @param mod model object
 * @param dat data object
 * @param opt options object
 * @return null
 */
double em_2_steps(model *mod, data *dat, options *opt)
{
	double L1;
	L1 = em_step(opt,dat,mod);
	mod->logL = L1;

	switch(mod->current_n) {
		case 0:
			COPY_3JAGGED_ARRAY(mod->init_pKLM,mod->pKLM,dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained)
				COPY_2ARRAY(mod->init_etaik, mod->etaik, mod->K);
			else
				COPY_1ARRAY(mod->init_etak, mod->etak, mod->K);
			break;
		case 1:
			COPY_3JAGGED_ARRAY(mod->iter1_pKLM,mod->pKLM,dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained)
				COPY_2ARRAY(mod->iter1_etaik, mod->etaik, mod->K);
			else
				COPY_1ARRAY(mod->iter1_etak, mod->etak, mod->K);
			break;
		case 2:
			COPY_3JAGGED_ARRAY(mod->iter2_pKLM, mod->pKLM, dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained)
				COPY_2ARRAY(mod->iter2_etaik, mod->etaik, mod->K);
			else
				COPY_1ARRAY(mod->iter2_etak, mod->etak, mod->K);
	}

	L1 = em_step(opt,dat,mod);
	mod->logL = L1;
	mod->current_n++;

	switch(mod->current_n){
		case 0:
			COPY_3JAGGED_ARRAY(mod->init_pKLM, mod->pKLM, dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained)
				COPY_2ARRAY(mod->init_etaik, mod->etaik, mod->K);
			else
				COPY_1ARRAY(mod->init_etak,mod->etak,mod->K);
			break;
		case 1:
			COPY_3JAGGED_ARRAY(mod->iter1_pKLM, mod->pKLM, dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained)
				COPY_2ARRAY(mod->iter1_etaik, mod->etaik, mod->K);
			else
				COPY_1ARRAY(mod->iter1_etak, mod->etak, mod->K);
			break;
		case 2:
			COPY_3JAGGED_ARRAY(mod->iter2_pKLM, mod->pKLM, dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained)
				COPY_2ARRAY(mod->iter2_etaik, mod->etaik, mod->K);
			else
				COPY_1ARRAY(mod->iter2_etak, mod->etak, mod->K);
	}
	mod->current_n++;
	return L1;
} /* em_2_steps */

/**
 * Print current parameter estimates.
 *
 * @param mod model object
 * @param dat data object
 * @param opt options object
 * @return null
 */
void print_param(options *opt, data *dat, model *mod, int which)
{
	int i, k, l, m;


	if (opt->admixture && !opt->eta_constrained) {
		for (i = 0; i < dat->I; i++) {
			fprintf(stderr, "etaik[%d]:", i);
			for (k = 0; k < mod->K; k++) {
				if (!which)
					fprintf(stderr, " %f", mod->etaik[i][k]);
				else if (which == 1)
					fprintf(stderr, " %f", mod->init_etaik[i][k]);
				else if (which == 2)
					fprintf(stderr, " %f", mod->iter1_etaik[i][k]);
				else if (which == 3)
					fprintf(stderr, " %f", mod->iter2_etaik[i][k]);
			}
			if (!((i+1) % 4))
				fprintf(stderr, "\n");
			else
				fprintf(stderr, " ");
		}
	} else {
		fprintf(stderr, "etaik:");
		for (k = 0; k < mod->K; k++)
			fprintf(stderr, " %f", mod->etak[k]);
		fprintf(stderr, "\n");
	}

	for (l = 0; l < dat->L; l++)
		for (k = 0; k < mod->K; k++) {
			fprintf(stderr, "Population %d, locus %d:", k, l);
			for (m = 0; m < dat->uniquealleles[l]; m++) {
				if (!which)
					fprintf(stderr, " %f", mod->pKLM[k][l][m]);
				else if (which == 1)
					fprintf(stderr, " %f", mod->init_pKLM[k][l][m]);
				else if (which == 2)
					fprintf(stderr, " %f", mod->iter1_pKLM[k][l][m]);
				else if (which == 3)
					fprintf(stderr, " %f", mod->iter2_pKLM[k][l][m]);
			}
			fprintf(stderr, "\n");
		}
} /* print_param */
