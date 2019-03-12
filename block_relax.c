/**
 * @file block_relax.c
 * @author Arun Sethuraman
 * @date Sat Jan 19 11:10:01 CST 2013
 * 
 * This file contains functions that will perform block relaxation updates
 * on mixing proportions and allele frequencies. This code does not work 
 * because it uses EM updates of the blocks; it should use sequential quadratic
 * programming to update the blocks.
 */

#include "multiclust.h"
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN /* return on memory allocation error */

/* block update */
double block_step(data *dat, model *mod, options *opt, double loglik);
double block_update_eta(options *opt, data *dat, model *mod, double prev_loglik);
double block_update_pkla(options *opt, data *dat, model *mod, double prev_loglik);

/* reset parameters after iteration */
void set_params_br(model *mod, data *dat, options *opt);
void retract_params(model *mod, data *dat, options *opt);

/**
 * Block relaxation (BR) algorithm.  This function attempts a full round of
 * block updates starting from the block after the last block that failed to
 * increase the likelihood.  If block update fails, model::em_done is set and
 * this algorithm aborts and returns the current log likelihood.
 *
 * @param dat data object
 * @param mod mod object
 * @param opt options object
 * @return log likelihood
 */
double block_update_params(options *opt, data *dat, model *mod)
{
	double loglik = 0;	/* log likelihood at next step */
//	int eta_updated = 1;	/* indicate if eta updated first */

	/* update eta or pkla block */
//	if (mod->em_done && (mod->current_l || mod->current_k)) {
		loglik = block_update_pkla(opt, dat, mod, mod->logL);
//		eta_updated = 0;
//	} else
//		loglik = block_update_eta(opt, dat, mod, mod->logL);

	if (isnan(loglik)) {	/* DEBUGGING */
		message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			CUSTOM_ERROR, "nan");
		exit(0);
	}

	/* resorted to EM update so iteration is done: return immediately */
	if (mod->em_done)
		return loglik;

	/* update other block */
//	if (eta_updated)
//		loglik = block_update_pkla(opt, dat, mod, loglik);
//	else
		loglik = block_update_eta(opt, dat, mod, loglik);

	if (isnan(loglik)) {	/* DEBUGGING */
		message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			CUSTOM_ERROR, "nan");
		exit(0);
	}

	return loglik;
} /* End of block_update_params */


/**
 * Block update of mixing proportions.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @param prev_loglik current log likelihood
 * @return current log likelhood
 */
double block_update_eta(options *opt, data *dat, model *mod, double prev_loglik)
{
	int debug = 0;
	int i;			/* current individual */
	double loglik = 0;

	mod->etaupdate = 1;	/* used by downstream functions */
	mod->em_done = 0;	/* reset; attempt accelerated updates */

	if (opt->admixture && !opt->eta_constrained) {
		i = 0;
		do {
			/* initialize finite differences for quasi-Newton update */
			initialize_etaiks(mod);

			/* block relax step on indiv. i: update model::etaik */
			loglik = block_step(dat, mod, opt, prev_loglik);

			if (isnan(loglik)) {	/* DEBUGGING */
				message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, CUSTOM_ERROR, "nan");
				exit(0);
			}

			/* block update produced smaller log likelihood */
			if (loglikelihood_decrease(opt, prev_loglik, loglik)) {
				if (debug) {
					fprintf(stderr, "Log likelihood "
						"decreased for etaik[%d] "
						"(%.15f < %.15f): ",
						mod->current_i, loglik,
						prev_loglik);
					for (int k = 0; k < mod->K; k++)
						fprintf(stderr, " %f",	
							mod->etaik[mod->current_i][k]);
					fprintf(stderr, " -> ");
				}

				/* retract to etaik update that failed */
				retract_params(mod, dat, opt);

				if (debug) {
					for (int k = 0; k < mod->K; k++)
						fprintf(stderr, " %f",
							mod->etaik[mod->current_i][k]);
					fprintf(stderr, "\n");
				}

				/* use EM to update all parameters */
				loglik = em_e_step(opt, dat, mod); 

				if (loglik < prev_loglik) {	/* DEBUGGING */
					message(stderr, __FILE__, __func__,
						__LINE__, ERROR_MSG, CUSTOM_ERROR,
						"ascent property violated: "
						"%f < %f\n", loglik, prev_loglik);
					exit(0);
				}

				/* update last iterate from EM update */
				set_params_em(mod, dat, opt);
				mod->em_done = 1;
				mod->current_i = (mod->current_i + 1) % dat->I;

				/* iteration done */
				return loglik;

			/* avoid nonsubstantive drop in log likelihood */
			} else if (prev_loglik > loglik) {

				/* reset log likelihood back to what it was */
				loglik = prev_loglik;

				/* also reset parameters to what they were */
				retract_params(mod, dat, opt);
			}

			set_params_br(mod, dat, opt);
			prev_loglik = loglik;

			mod->current_i = (mod->current_i + 1) % dat->I;
		} while (++i < dat->I);
exit(0);
	
	/* This procedure runs BR on the mixture model etaks or mixing proportions */
	} else  {

		/* Initialize all the etaks */
		initialize_etaks(mod);

		/* Run BR algorithm */
		loglik = block_step( dat, mod, opt, prev_loglik);

		if (isnan(loglik)) {	/* DEBUGGING */
			message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, CUSTOM_ERROR, "nan");
			exit(0);
		}

		if (loglikelihood_decrease(opt, prev_loglik, loglik)) {
			if (debug)
				fprintf(stderr, "Log likelihood decreased for "
					"eta (%.15f < %.15f)!\n", loglik,
					prev_loglik);
			/* Resorting to EM if the likelihood isn't improving */
			retract_params(mod, dat, opt);
			loglik = em_e_step(opt, dat, mod);

/* [KSD] [TODO] [BUG] [WORKING] Something wrong here! */
			if (loglik - prev_loglik < -1e-7) {	/* DEBUGGING */
				message(stderr, __FILE__, __func__, __LINE__,
					ERROR_MSG, CUSTOM_ERROR, "ascent "
					"property violated: %f < %f (%e)\n", 
					loglik, prev_loglik,
					loglik - prev_loglik);
				exit(0);
			}

			set_params_em(mod,dat,opt); /* Set all updated parameters post EM */
			mod->em_done = 1;
			return loglik;
		} else if (prev_loglik > loglik) {
			loglik = prev_loglik;
			retract_params(mod, dat, opt);
		}

		set_params_br(mod,dat,opt);
		prev_loglik = loglik;
	}

	return loglik;
} /* block_update_eta */


/**
 * Block update of pkla.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @param prev_loglik current log likelihood
 * @return current log likelihood
 */
double block_update_pkla(options *opt, data *dat, model *mod, double prev_loglik)
{
	int debug = 0;
	double loglik = 0;	/* log likelihood at next step */
	int l = 0;		/* l indexes current locus */
	int stop = dat->L * mod->K;


	mod->etaupdate = 0;	/* [KSD] Used by downstream functions */
	mod->em_done = 0;

	if (debug>1)
		fprintf(stderr, "block_update_pkla(%f)\n", prev_loglik);

	do {
		if (debug>1)
			fprintf(stderr, "block_update_pkla(%d,%d)\n",
				mod->current_k, mod->current_l);

		/* Initialize subpopulation allele frequencies at locus l*/
		initialize_pklas(mod, dat);

		/* block update for this subpopulation, locus combo */
		loglik = block_step(dat, mod, opt, prev_loglik);

		if (isnan(loglik)) {	/* DEBUGGING */
			message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
				CUSTOM_ERROR, "nan");
			exit(0);
		}
			
		/* Skipping into EM if likelihood is not improving*/
		if (loglikelihood_decrease(opt, prev_loglik, loglik)) {
			if (debug) {
				fprintf(stderr, "Log likelihood decreased for "
					"p[%d][%d]a (%.15f < %.15f):",
					mod->current_k, mod->current_l, loglik,
					prev_loglik);
				for (int m = (dat->L_alleles
					&& dat->L_alleles[mod->current_l][0]
					== MISSING ? 1 : 0);
					m < dat->uniquealleles[mod->current_l];
					m++)
					fprintf(stderr, " %f",
						mod->pKLM[mod->current_k][
						mod->current_l][m]);
				fprintf(stderr, " -> ");
			}

			retract_params(mod, dat, opt);

			/* Call EM step to update all parameters */
			loglik = em_e_step(opt, dat, mod);

			if (debug) {
				for (int m = (dat->L_alleles
					&& dat->L_alleles[mod->current_l][0]
					== MISSING ? 1: 0);
						 m < dat->uniquealleles[mod->current_l]; m++)
					fprintf(stderr, " %f",
						mod->pKLM[
						mod->current_k][
						mod->current_l][m]);
				fprintf(stderr, "\n");
			}

			/* Update all parameters post EM */
			set_params_em(mod, dat, opt);
			mod->em_done = 1;
			mod->current_l = (mod->current_l + 1) % dat->L;
			if (mod->current_l == 0)
				mod->current_k = (mod->current_k + 1) % mod->K;
			return loglik;
		} else if (prev_loglik > loglik) {
			loglik = prev_loglik;
			retract_params(mod, dat, opt);
		}
		
		/* Set parameters from BR estimates and continue */
		set_params_br(mod, dat, opt);
		prev_loglik = loglik;
		mod->current_l = (mod->current_l + 1) % dat->L;
		if (mod->current_l == 0)
			mod->current_k = (mod->current_k + 1) % mod->K;
	} while (++l < stop);

	return loglik;
} /* End of block_update_pkla */


/**
 * Block update of a subset of parameters. This function is called by
 * block_update_params() to update mixing proportions, etaiks, of one
 * individual, i, or to update allele frequencies, pklas, at a locus l.
 *
 * @param dat data object
 * @param opt option object
 * @param mod mod object
 * @return double log likelihood
 */
double block_step(data *dat, model *mod, options *opt, double loglik)
{
	int debug = 1;
	double numer = 0.0;
	double denom = 0.0;
	double S;
	double sum = 0.0;
	int below_0 = 0, above_1 = 0;
	int i, l = 0, k, k1, m, m_start = 0;
	double lterm;

	i = 0;
	if (mod->etaupdate == 1) {
	
		//i = mod->current_i;
		for (k = 0; k < mod->K; k++) {	
			numer += mod->U[k] * mod->U[k];
			denom += (mod->V[k] - mod->U[k]) * (mod->V[k] - mod->U[k]);
		}

		S = (-1)*sqrt(numer/denom);

		if (isnan(S))
			return loglik;

		if (debug && fabs(S) > 100)
			fprintf(stderr, "S = %f\n", S);
	
		if (opt->admixture && !opt->eta_constrained) {
			i = mod->current_i;
			/* remove part of log likelihood due to this individual */
			for (l = 0; l < dat->L; l++) {
				m_start = dat->L_alleles
					&& dat->L_alleles[l][0] == MISSING
					? 1 : 0;
				for (m = m_start; m < dat->uniquealleles[l]; m++) {
					if (!dat->ILM[i][l][m])
						continue;
					lterm = 0;
					for (k = 0; k < mod->K; k++)
						lterm += mod->etaik[i][k] * mod->pKLM[k][l][m];
					loglik -= dat->ILM[i][l][m] * log(lterm);
				}
			}

			/* Calculating new admixture proportions using SqS3 update */
			for (k = 0; k < mod->K; k++) {
				mod->etaik[i][k] = mod->etaik[i][k]
					- 2*S*mod->U[k]
					+ S*S*(mod->V[k]- mod->U[k]); 

				/* [KSD] Does the algorithm say this?
 				[AS] No - but I added it to prevent it from going -ve.
				I suppose it is checked by the projection too so might
				not be needed. 				
				if (mod->etaik[i][k] < 0.0)
					mod->etaik[i][k] = 0.0;
				*/

				/* check on constraint */
				sum += mod->etaik[i][k];
				if (mod->etaik[i][k] < 0)
					below_0 = 1;
				else if (mod->etaik[i][k] > 1)
					above_1 = 1;
			}

		} else {
			for(k = 0; k < mod->K; k++){
				mod->etak[k] = mod->etak[k]
					- 2*S*mod->U[k]
					+ S*S*(mod->V[k] - mod->U[k]);
				/* [KSD] Does the algorithm say this?
 				[AS] - see above.		
				if(mod->etak[k] < 0.0)
					mod->etak[k] = 0.0;
				*/

				/* check on constraint */
				sum += mod->etak[k];
				if (mod->etak[k] < 0)
					below_0 = 1;
				else if (mod->etak[k] > 1)
					above_1 = 1;
			}
		}

	} else {
		l = mod->current_l;
		k = mod->current_k;

		/* [KSD] [TODO] How should we deal with missing? */
		m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
                for (m = m_start; m < dat->uniquealleles[l]; m++) {
			numer += mod->U[m] * mod->U[m];
			denom += mod->V[m] - mod->U[m] * (mod->V[m] - mod->U[m]);
		}

		S = (-1)*sqrt(numer/denom);
if (debug)
fprintf(stderr, "%d, %d (pre): %f", l, k, S);
		if (isnan(S))
			return loglik;

		/* remove part of log likelihood affected by this block */
		for (i = 0; i < dat->I; i++)
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				if (!dat->ILM[i][l][m])
					continue;
				lterm = 0;
				for (k1 = 0; k1 < mod->K; k1++)
					lterm += (opt->admixture
						&& !opt->eta_constrained
						? mod->etaik[i][k1]
						: mod->etak[k1]) * mod->pKLM[k1][l][m];
				loglik -= dat->ILM[i][l][m] * log(lterm);
			}

               	for (m = m_start; m < dat->uniquealleles[l]; m++) {
			mod->pKLM[k][l][m] = mod->pKLM[k][l][m] - 2*S*mod->U[m]
				+ S*S*(mod->V[m] - mod->U[m]);
			/* [KSD] Does the algorithm say this??
			if (mod->pKLM[k][l][m] < 0.0)
				mod->pKLM[k][l][m] = 0;
			*/
if (debug)
fprintf(stderr, " %f", mod->pKLM[k][l][m]);
			/* check on constraint */
			sum += mod->pKLM[k][l][m];
			if (mod->pKLM[k][l][m] < 0)
				below_0 = 1;
			else if (mod->pKLM[k][l][m] > 1)
				above_1 = 1;
		}
if (debug)
fprintf(stderr, "\n");
	
	}

	/* If admixture or allele proportions don't sum to 1, project onto simplex */
	/* [KSD] pretty stringent requirement; is this costing us? */
	if (fabs(sum - 1.0) > 2*DBL_EPSILON || above_1 || below_0) {
		if (debug) {
			fprintf(stderr, "simplex_project(%e, %d, %d):",
				sum - 1.0, above_1, below_0);
			if (mod->etaupdate && opt->admixture && !opt->eta_constrained) {
				for (k = 0; k < mod->K; k++)
					fprintf(stderr, " %f", mod->etaik[i][k]);
			} else if (mod->etaupdate) {
				for (k = 0; k < mod->K; k++)
					fprintf(stderr, " %f", mod->etak[k]);
			} else {
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					fprintf(stderr, " %f", mod->pKLM[k][l][m]);
			}
			fprintf(stderr, "\n");
		}
		if (mod->etaupdate)
			simplex_project_eta(mod, opt, mod->current_i);
		else
			simplex_project_pklm(mod, dat, opt, mod->current_k, mod->current_l);
	}

fprintf(stderr, "%d, %d (pos):", l, k);
for (m = m_start; m < dat->uniquealleles[l]; m++)
fprintf(stderr, " %f", mod->pKLM[k][l][m]);
fprintf(stderr, "\n");

	/* add updated part of log likelihood due to this block */
	if (mod->etaupdate && opt->admixture && !opt->eta_constrained) {
		for (l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING ? 1 : 0;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {

				if (!dat->ILM[i][l][m])
					continue;

				lterm = 0;
				for (k = 0; k < mod->K; k++)
					lterm += mod->etaik[i][k] * mod->pKLM[k][l][m];
				loglik += dat->ILM[i][l][m] * log(lterm);

				if (isnan(loglik)) {	/* DEBUGGING */
					message(stderr, __FILE__, __func__,
						__LINE__, ERROR_MSG,
						CUSTOM_ERROR, "l = %d, m = %d: "
						"lterm = %f, loglik = %f\n", l,
						m, lterm, loglik);
					exit(0);
				}
			}
		}

		return loglik;
	} else if (!mod->etaupdate) {
		for (i = 0; i < dat->I; i++) {
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING ? 1 : 0;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				if (!dat->ILM[i][l][m])
					continue;
				lterm = 0;
				for (k1 = 0; k1 < mod->K; k1++)
					lterm += (opt->admixture && !opt->eta_constrained ? 
						mod->etaik[i][k1] : mod->etak[k1]) * mod->pKLM[k1][l][m];
					//lterm += mod->etaik[i][k1] * mod->pKLM[k1][l][m];		
				loglik += dat->ILM[i][l][m] * log(lterm);

				if (isnan(loglik)) {	/* DEBUGGING */
					message(stderr, __FILE__, __func__,
						__LINE__, ERROR_MSG,
						CUSTOM_ERROR, "i = %d, m = %d: "
						"lterm = %f, loglik = %f\n", i,
						m, lterm, loglik);
					exit(0);
				}
			}
		}

		return loglik;
	}

	/* [KSD] We are here if mod->etaupdate && opt->eta_constrained, in
	 * which case we have to compute the whole log likelihood again. */
	return log_likelihood(opt, dat, mod, 0);
} /* block_step */


/**
 * Compute last two finite differences in pkla block.
 *
 * @param mod model object
 * @param dat data object
 * @return none
 */
void initialize_pklas(model *mod, data *dat)
{
	int k, l, m, m_start;

	l = mod->current_l;
	k = mod->current_k;
	m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;

	for (m = m_start; m < dat->uniquealleles[l] ; m++) {
		/* init_ -> iter1_ -> iter2_ -> init_ -> iter1_ -> ... */
		switch (mod->current_n) {
			case 0:	
				mod->U[m] = mod->iter1_pKLM[k][l][m]
					- mod->init_pKLM[k][l][m];
				mod->V[m] = mod->iter2_pKLM[k][l][m]
					- mod->iter1_pKLM[k][l][m];
				break;
			case 1:
				mod->U[m] = mod->iter2_pKLM[k][l][m]
					- mod->iter1_pKLM[k][l][m];
				mod->V[m] = mod->init_pKLM[k][l][m]
					- mod->iter2_pKLM[k][l][m];
				break;
			case 2:
				mod->U[m] = mod->init_pKLM[k][l][m]
					- mod->iter2_pKLM[k][l][m];
				mod->V[m] = mod->iter1_pKLM[k][l][m]
					- mod->init_pKLM[k][l][m];
				break;
		}
	}				
}/* End of initialize_pklas */


/**
 * Compute last two finite differences in etaik block.
 *
 * @param mod model object
 * @return none
 */
void initialize_etaiks(model *mod)
{
	int k;
	int i;

	/* current individual */
	i = mod->current_i;

	for (k = 0; k < mod->K ; k++) {
		switch (mod->current_n) {
			case 0:		
				mod->U[k] = mod->iter1_etaik[i][k]
					- mod->init_etaik[i][k];
				mod->V[k] = mod->iter2_etaik[i][k]
					- mod->iter1_etaik[i][k];
				break;	/* [KSD] was critical bug?? */
			case 1:
				mod->U[k] = mod->iter2_etaik[i][k]
					- mod->iter1_etaik[i][k];
				mod->V[k] = mod->init_etaik[i][k]
					- mod->iter2_etaik[i][k];
				break;
			case 2:
				mod->U[k] = mod->init_etaik[i][k]
					- mod->iter2_etaik[i][k];
				mod->V[k] = mod->iter1_etaik[i][k]
					- mod->init_etaik[i][k];
				break;
		}				
	}
} /* End of initialize_etaiks function */


/**
 * Compute last two finite differences in the etak block.
 *
 * @param mod model object
 * @return none
 */
void initialize_etaks(model *mod)
{
	int k;

	for (k = 0; k < mod->K; k++) {
		switch (mod->current_n) {
			case 0:		
				mod->U[k] = mod->iter1_etak[k]
					- mod->init_etak[k];
				mod->V[k] = mod->iter2_etak[k]
					- mod->iter1_etak[k];
				break; /* [KSD] was BUG */
			case 1:
				mod->U[k] = mod->iter2_etak[k]
					- mod->iter1_etak[k];
				mod->V[k] = mod->init_etak[k]
					- mod->iter2_etak[k];
				break;	/* [KSD] was BUG */
			case 2:
				mod->U[k] = mod->init_etak[k]
					- mod->iter2_etak[k];
				mod->V[k] = mod->iter1_etak[k]
					- mod->init_etak[k];
				break;
		}
	}
} /* End of initialize_etaks function */


/**
 * Update iterates after a successful accelerated block update.  This function
 * is updating the next iterate block-by-block as successful acceleration steps
 * are completed.  If a accelerated block update fails, then the the algorithm
 * falls back on an EM iteration, which will wholesale replace the next iterate.
 * See set_params_em() for the function that does this wholesale replace.
 *
 * @param mod model object
 * @param opt options object
 * @param dat data object
 * @return int
 */
void set_params_br(model *mod, data __attribute__((unused)) *dat, options *opt)
{
	int i, k, l, m;
	int m_start;

	if (opt->admixture && !opt->eta_constrained && mod->etaupdate == 1) {
		i = mod->current_i;
		switch (mod->current_n) {
			case 0:
				for (k = 0; k < mod->K; k++)
					mod->init_etaik[i][k] = mod->etaik[i][k];
				break;
			case 1:
				for (k = 0; k < mod->K; k++)
					mod->iter1_etaik[i][k] = mod->etaik[i][k];
				break;
			case 2:
				for (k = 0; k < mod->K; k++)
					mod->iter2_etaik[i][k] = mod->etaik[i][k];
				break;
		}
	} else if (mod->etaupdate == 1) {
		switch (mod->current_n) {
			case 0:
				for (k = 0; k < mod->K; k++)
					mod->init_etak[k] = mod->etak[k];
				break;
			case 1:
				for (k = 0; k < mod->K; k++)
					mod->iter1_etak[k] = mod->etak[k];
				break;
			case 2:
				for (k = 0; k < mod->K; k++)
					mod->iter2_etak[k] = mod->etak[k];
				break;
		}
	} else {
		l = mod->current_l;
		k = mod->current_k;

		/* [KSD] [TODO] how to deal with missing alleles */
		m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING
			? 1 : 0;

		switch (mod->current_n) {
			case 0:
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					mod->init_pKLM[k][l][m] = mod->pKLM[k][l][m];
				break;
			case 1:
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					mod->iter1_pKLM[k][l][m] = mod->pKLM[k][l][m];
				break;
			case 2:
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					mod->iter2_pKLM[k][l][m] = mod->pKLM[k][l][m];
				break;
		}
	}
} /* End of set_params_br function */


/**
 * Update block accerated iterates after EM step.  This function updates
 * everything if we are using an EM during the block relaxation algorithm
 *
 * @param mod model object
 * @param dat data object
 * @param opt options object
 * @return state
 */
int set_params_em(model *mod, data *dat, options *opt)
{
	switch (mod->current_n) {
		case 0:
			COPY_3JAGGED_ARRAY(mod->init_pKLM, mod->pKLM, dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained) {
				COPY_2ARRAY(mod->init_etaik, mod->etaik, mod->K);
			} else
				COPY_1ARRAY(mod->init_etak, mod->etak, mod->K);
			break;
		case 1:
			COPY_3JAGGED_ARRAY(mod->iter1_pKLM, mod->pKLM, dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained) {
				COPY_2ARRAY(mod->iter1_etaik, mod->etaik, mod->K);
			} else
				COPY_1ARRAY(mod->iter1_etak, mod->etak, mod->K);
			break;
		case 2:
			COPY_3JAGGED_ARRAY(mod->iter2_pKLM, mod->pKLM, dat->uniquealleles);
			if (opt->admixture && !opt->eta_constrained) {
				COPY_2ARRAY(mod->iter2_etaik, mod->etaik, mod->K);
			} else
				COPY_1ARRAY(mod->iter2_etak, mod->etak, mod->K);
			break;
	}

	return mod->current_n;
} /* End of set_params_em function */


/**
 * Resets the parameters after failed accelerated update.  If the accelerated
 * update fails, there next iterate may be partially updated.  However, the
 * last block updated has not improved the likelihood and needs to be retracted
 * and replaced with the last estimates.  The last estimates are stored in
 * model::init_*, model::iter1_*, or model::iter2_* variables, depending on the
 * value of model::current_n.
 *
 * @param mod model object
 * @param opt options object
 * @return null
 */
void retract_params(model *mod, data __attribute__((unused)) *dat, options *opt)
{
	int i, k, l, m, m_start;

	/* [KSD] [TODO] how to handle missing data? */
	
	if (mod->etaupdate == 1 && opt->admixture && !opt->eta_constrained) {
		i = mod->current_i;
		switch (mod->current_n) {
			case 0:
				for(k = 0; k < mod->K; k++)
					mod->etaik[i][k]
						= mod->iter2_etaik[i][k];
				break;
			case 1:
				for(k = 0; k < mod->K; k++)
					mod->etaik[i][k]
						= mod->init_etaik[i][k];
				break;
			case 2:
				for(k = 0; k < mod->K; k++)
					mod->etaik[i][k]
						= mod->iter1_etaik[i][k];
				break;
		}
	} else if (mod->etaupdate == 1) {
		/* reverse attempted update of etak */
		switch (mod->current_n) {
			case 0:
				for (k = 0; k < mod->K; k++)
					mod->etak[k] = mod->iter2_etak[k];
				break;
			case 1:
				for (k = 0; k < mod->K; k++)
					mod->etak[k] = mod->init_etak[k];
				break;
			case 2:
				for (k = 0; k < mod->K; k++)
					mod->etak[k] = mod->iter1_etak[k];
				break;
		}
	} else {
		/* return to previous update of unsuccessful pKLM update */
		k = mod->current_k;
		l = mod->current_l;
		m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING
			? 1 : 0;
		switch(mod->current_n) {
			case 0:
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					mod->pKLM[k][l][m] = mod->iter2_pKLM[k][l][m];
				break;
			case 1:
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					mod->pKLM[k][l][m] = mod->init_pKLM[k][l][m];
				break;
			case 2:
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					mod->pKLM[k][l][m] = mod->iter1_pKLM[k][l][m];
				break;
		}
	}
} /* End of retract_params function */			
