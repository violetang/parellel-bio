/**
 * @file multiclust.c
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Wed Dec  5 09:20:48 CST 2012
 *
 * User interface to multiclust.
 *
 * TODOS
 * [TODO] write the results of parametric bootstrap to a permanent output file
 * [TODO] run matrix of hypothesis tests and control FDR (as per Maitra & 
 *        Melnykov) to auto-estimate K
 * [TODO] add hypothesis test for H0: mixture vs. HA: admixture
 * [TODO] add code to run variable number of initializations, perhaps until the
 *        current best log likelihood has been observed X times
 */

#include "mpi.h"
#include <time.h>
#include "cline.h"	/* command line */
#include "multiclust.h"
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN

/* create structures for options, data, and model */
int make_options(options **opt);
int make_data(data **dat);
int make_model(model **);

/* set/parse options and model; verify integrity */
int parse_options(options *opt, data *dat, int argc, const char **argv);
int allocate_model(options *opt, model *mod, data *dat);
int synchronize(options *opt, data *dat, model *mod);

/* do the data fitting or bootstrap */
int estimate_model(options *opt, data *dat, model *mod, int bootstrap);
int maximize_likelihood(options *opt, data *dat, model *mod, int bootstrap);
int run_bootstrap(options *opt, data *dat, model *mod);

/* cleanup all allocated memory */
void free_options(options *opt);
void free_data(data *dat);
void free_model(model *mod, options *opt);
void free_model_mles(model *mod);
void free_model_data(model *mod, options *opt);

const char *accel_method_abbreviations[NUM_ACCELERATION_METHODS] = {
        "",
	"S1",
	"S2",
	"S3",
	"Q1"
};
const char *accel_method_names[NUM_ACCELERATION_METHODS] = {
	"No acceration",
	"SQUAREM version 1",
	"SQUAREM version 2",
	"SQUAREM version 3",
	"Quasi Newton (q=1)"
};




int main(int argc, const char **argv)
{
	options *opt = NULL;	/* run options */
	data *dat = NULL;	/* genetic data */
	model *mod = NULL;	/* model parameters */
	int err = NO_ERROR;	/* error code */

	/* make various structures to store run information */
	if ((err = make_options(&opt)))
		goto FREE_AND_EXIT;

	/* set up data structure */
	if ((err = make_data(&dat)))
		goto FREE_AND_EXIT;

	/* set up model structure */
	if ((err = make_model(&mod)))
		goto FREE_AND_EXIT;

	/* parse command-line options */
	if ((err = parse_options(opt, dat, argc, argv)))
		goto FREE_AND_EXIT;


	/*
	* MPI environment set up
	*/

	int my_rank;
	//int p = opt-> process;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &opt->process);
	float dest = 0;
	int tag = 0;
	int source;
	MPI_Status status;

	/* read data */
	if ((err = read_file(opt, dat)))
		goto FREE_AND_EXIT;

	/* finalize settings that refer to data; check settings */
	if ((err = synchronize(opt, dat, mod)))
		goto FREE_AND_EXIT;

	/* estimate the model(s) using the observed data */
	if ((err = estimate_model(opt, dat, mod, 0)))
		goto FREE_AND_EXIT;

	/* optionally run a bootstrap */
	if (opt->n_bootstrap) {
		if ((err = run_bootstrap(opt, dat, mod)))
			goto FREE_AND_EXIT;

		/* TODO write this result to some file! */
		fprintf(stdout, "p-value to reject H0: K=%d is %f\n", 
			mod->null_K, mod->pvalue);
	}

FREE_AND_EXIT:

	free_model(mod,opt);
	free_options(opt);
	free_data(dat);

	MPI_Finalize();

	return err;

} /* main */


/**
 * Estimate model(s) using maximum likelihood.  For each model, the function
 * allocates the model (parameters), maximizes the likelihood (probably using
 * multiple initializations), extracts relevant statistics (mles, log
 * likelihood), and frees the model.  Currently, the code either tests H0
 * against HA (and aborts if the HA log likelihood does not exceed the H0 log
 * likelihood) or fits a series of models from K=1 to K=_options::max_K.
 * Importantly, when hypothesis testing, it fits H0 //before// HA.  One could
 * modify this code to select other models for fitting.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @param bootstrap indicate if bootstrap run
 * @return error status
 */
int estimate_model(options *opt, data *dat, model *mod, int bootstrap)
{
	int total_iter = 0;
	int err = NO_ERROR;

	/* initialize: if bootstrap, comparing H0: K = mod->null_K vs.
	 * Ha: K = mod->alt_K; otherwise running for 1, 2, ..., opt->max_K
	 * so log likelihood will increase as progress through models */
	mod->max_logL = -Inf;
	mod->K = opt->n_bootstrap ? mod->null_K : opt->min_K;

	//if(opt->block_relax == 1)
	//	mod->K = 2;

	do { //when not parallel, only one processor, so need loop to run agian with different start random possibility
		 //if parallel, do we still need use loop?
		
		
		/* Setting the largest size of the U and V arrays to be the larger of K and max_M */
		if (dat->max_M < mod->K)
			dat->max_M = mod->K;
		
		/* knowing K, can allocate space for parameters */
		if ((err = allocate_model(opt, mod, dat)))
			return err;

		/* maximize the likelihood; involves multiple initializations */
		/* also stores mles in mod->mle_* parameters */
		if ((err = maximize_likelihood(opt, dat, mod, bootstrap)))
			return err;

		total_iter += mod->n_total_iter;
		//here you would have mod->max_logL for each processor
		//From each processor, send that max_logL to the head node

		/* possibly store maximum likelihood under H0 */
		if (opt->n_bootstrap && mod->K == mod->null_K)
			mod->max_logL_H0 = mod->max_logL;

		/* was BUG: popq_admix() and indivq_admix() written here for
		 * last iteration, not best iteration */

		//after get the max possibility, the processor will clean all the data.
		//so before free all the data, need save each processor's max possibility here
		//then compare all processors' max ppssibility, only save the biggiest one. 

		/* free parameters */
		free_model_data(mod,opt);

		/* resetting some model params */
		mod->current_i = 0;
		mod->current_l = 0;
		mod->current_k = 0;
		mod->current_n = 0;
		mod->logL = 0.0;
		mod->converged = 0;
		mod->etaupdate = 0;
		mod->em_done = 0;
		/* choose next K */
		if (opt->n_bootstrap && mod->K == mod->null_K)
			mod->K = mod->alt_K;
		else if (!opt->n_bootstrap && mod->K < opt->max_K)
			mod->K++;
		else
			break;
	} while (1);

	if (opt->n_init > 1)
		fprintf(stderr, "Total number of iterations: %d\n", total_iter);

	/* if running bootstrap, record test statistic */
	if (opt->n_bootstrap) {
		double diff = mod->max_logL - mod->max_logL_H0;

		if (diff <= 0)	/* == OK, but prob. convergence error */
			return message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, INTERNAL_ERROR, "Null hypothesis "
				"likelihood exceeds alternative hypothesis "
				"likelihood.  Try increasing number of "
				"initializations (command-line option -n)\n");

		if (!bootstrap)
			mod->ts_obs = diff;
		else 
			mod->ts_bs = diff;
	}

	return err;
} /* estimate_model */


/**
 * Maximize the log likelihood.  Given a particular model (admixture/mixture
 * and _model::K), initialize the model _options::n_init initializations times,
 * maximize the likelihood, and record the result if a better solution is
 * found.  The best log likelihood is stored in _model::max_logL.  The mle
 * parameter estimates are stored in _model::mle_pKLM and _model::mle_etak (or
 * _model::mle_etaik), but only if the observed data is being fit to H0.  In
 * addition, if fitting the observed data, the best fitting solution is written
 * to file.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @param bootstrap is this a bootstrap run
 * @return error status
 */
int maximize_likelihood(options *opt, data *dat, model *mod, int bootstrap)
{
	int i;
	int err = NO_ERROR;
	int ntimes = 0;
	clock_t clk;

	for (i=0; opt->n_seconds || i<opt->n_init; i++) {

		/* initialize parameters */
		if ((err = initialize_model(opt, dat, mod)))
			return err;

		/* maximize likelihood */
		if (opt->accel_scheme) 
			accelerated_em(opt, dat, mod);
		else
			em(opt, dat, mod);

		mod->n_total_iter += mod->n_iter;

		/* better solution */
		if (mod->logL > mod->max_logL) {
			ntimes = 0;
			mod->max_logL = mod->logL;

			/* save mles if not bootstrap run and this is H0 */
			if (!bootstrap && opt->n_bootstrap && mod->K == mod->null_K) {
				COPY_3JAGGED_ARRAY(mod->mle_pKLM, mod->pKLM,
					dat->uniquealleles);
				if (!opt->admixture || opt->eta_constrained)
					COPY_1ARRAY(mod->mle_etak, mod->etak, mod->K);
				else
					COPY_2ARRAY(mod->mle_etaik, mod->etaik, mod->K);
			}

			/* write results to file if not bootstrap run */
			if (!bootstrap) {
				/* TODO [KSD]: overwriting potentially many 
				 * times for big data is bad */
				if (opt->admixture) {
					partition_admixture(dat, mod);
					write_file_detail(opt, dat, mod);
					popq_admix(opt, dat, mod);
					indivq_admix(opt, dat, mod);
				} else {
					partition_mixture(dat, mod);
					write_file_detail(opt, dat, mod);
					popq_mix(opt, dat, mod);
					indivq_mix(opt, dat, mod);
				}
			
			}
		} else if (stop_iterations(opt, mod, mod->max_logL))
			ntimes++;

		if (!bootstrap)
			fprintf(stdout, "K = %d, initialization = %d: %f "
				"in %3d iterations (%f; %d)\n", mod->K, i,
				mod->logL, mod->n_iter, (float) mod->max_logL,
				ntimes);
		
		if (mod->K != 1) {
			mod->current_i = 0;
			mod->current_l = 0;
			mod->current_k = 0;
			mod->current_n = 0;
			mod->logL = 0.0;
			mod->converged = 0;
			mod->etaupdate = 0;
			mod->em_done = 0;
		}

		/* global mle if only one cluster; no need for multiple
		 * initializations */
		if (mod->K == 1)
			break;

		if (opt->n_seconds) {
			clk = clock();
			if (clk < 0)
				return message(stderr, __FILE__, __func__,
					__LINE__, ERROR_MSG, INTERNAL_ERROR,
					"Timing information is not available "
					"on this computer, so -t option "
					"cannot be used.\n");
			else if (clk/CLOCKS_PER_SEC > opt->n_seconds)
				break;
		}
	}

	return err;
} /* End of maximize_likelihood(). */


/**
 * Parametric bootstrap.  For __options::n_bootstrap times, simulate bootstrap
 * dataset using parameter estimates stored in _model::mle_etak or
 * _model::mle_etaik and _model::mle_pKLM.  Fit the H0 and HA models using
 * exactly the same procedure used to fit the observed data to the same models.
 * Record the number of times the test statistic, stored in _model::ts_bs,
 * exceeds the observed test statistic, stored in _model::ts_obs.  Store the
 * resulting estimated p-value in _model::pvalue.  While running, this function
 * doubles the amount of memory required to store the data, but it clears the
 * memory when done.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int run_bootstrap(options *opt, data *dat, model *mod)
{
	int i;
	int ntime = 0;
	int err = NO_ERROR;

	for (i=0; i < opt->n_bootstrap; i++) {
		fprintf(stdout, "Bootstrap dataset %d (of %d):", i+1,
			opt->n_bootstrap);

		/* generate bootstrap dataset under H0 */
		if ((err = parametric_bootstrap(opt, dat, mod)))
			return err;
/* temporary : if you want to generate some simulation data to play with
write_data(opt, dat, 1);
*/

		/* fit models H0 and HA */
		if ((err = estimate_model(opt, dat, mod, 1)))
			return err;

		/* mod->max_logL is maximum log likelihood under HA */
		if (mod->ts_bs >+ mod->ts_obs)
			ntime++;
		fprintf(stdout, " test statistics bs=%f obs=%f (%f)\n",
			mod->ts_bs, mod->ts_obs, (double) ntime/(i+1));
	}

	mod->pvalue = ntime / opt->n_bootstrap;

	cleanup_parametric_bootstrap(dat);

	return err;
} /* run_bootstrap */


/**
 * Synchronize options, data, and model.  This function checks to make sure
 * that the user has made self-consistent choices.  The user cannot fit K
 * subpopulations when there are fewer than K observations.  The current
 * bootstrap hypothesis testing procedure reads -k <k> from the command line
 * and tests H0: K = <k> - 1 vs. HA: K = <k>.  This choice is checked and set
 * up here.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int synchronize(options *opt, data *dat, model *mod)
{
	opt->lower_bound = MIN(opt->lower_bound, 1.0 / dat->I / dat->ploidy - 0.5 / dat->I / dat->ploidy);
	if (dat->I < opt->max_K)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INVALID_USER_SETUP, "Maximum number of clusters (%d) "
			"(set with command-line argument -k) cannot exceed "
			"the number of individuals (%d)\n", opt->max_K, dat->I);
	if (opt->n_bootstrap && opt->max_K <= 1)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INVALID_USER_SETUP, "When bootstrapping, maximum K (%d) "
			"(set with command-line argument -k) must exceed 1.",
			opt->max_K);
	if (opt->n_bootstrap) {
		mod->null_K = opt->max_K - 1;
		mod->alt_K = opt->max_K;
	}

	if (opt->min_K > opt->max_K)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INVALID_USER_SETUP, "Minimum K (%d) must not exceed "
			"maximum K (%d).", opt->min_K, opt->max_K);

	return NO_ERROR;
} /* synchronize */


/**
 * Allocate options object and initialize.
 *
 * @param opt options object pointer reference
 * @return error status
 */
int make_options(options **opt)
{

	*opt = malloc(sizeof **opt);
	if (*opt == NULL)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, "options object");
	
	(*opt)->filename = NULL;
	(*opt)->R_format = 0;
	(*opt)->alleles_are_indices = 0;
	(*opt)->seed = 1234567;
	(*opt)->n_init = 50;

	(*opt)->process = 1 ; /*defalut process = 1*/

	/* convergence criterion set to match Lange's definition */
	(*opt)->max_iter = 0;	/*2000;*/
	(*opt)->rel_error = 0;	/*1e-6;*/
	(*opt)->abs_error = 1e-4;

	/* default tests K=6 */
	(*opt)->min_K = 6;
	(*opt)->max_K = 6;

	/* default initialization method: */
	(*opt)->initialization_method = RANDOM_CENTERS;	// TESTING;
	(*opt)->initialization_procedure = NOTHING;/*RAND_EM;*/
	(*opt)->n_rand_em_init = 50;
	(*opt)->lower_bound = 1e-8;
	(*opt)->path = "./";
	(*opt)->admixture = 0;
	(*opt)->eta_constrained = 0;
	(*opt)->n_bootstrap = 0;
	(*opt)->block_relax = 0;
	(*opt)->accel_scheme = 0;
	(*opt)->n_init_iter = 0;
	(*opt)->n_seconds = 0;
	(*opt)->fresh_em = 1;		/* Varadhan2008 */
	(*opt)->adjust_step = 0;//1000;	/* Varadhan2008: \infty */
	(*opt)->verbosity = SILENT;
	(*opt)->qfile = NULL;
	(*opt)->pfile = NULL;

	return NO_ERROR;
} /* make_options */


/**
 * Free options object.
 *
 * @param opt options object
 * @return void
 */
void free_options(options *opt)
{
	if (opt)
		free(opt);
	opt = NULL;
} /* free_options */


/**
 * Allocate data object and initialize.
 *
 * @param dat data object pointer reference
 * @return error status
 */
int make_data(data **dat)
{

	*dat = malloc(sizeof **dat);
	if (*dat == NULL)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, "data object");
	
	(*dat)->I = 0;
	(*dat)->L = 0;
	(*dat)->M = 0;
	(*dat)->ploidy = 2;
	(*dat)->IL = NULL;
	(*dat)->ila = NULL;
	(*dat)->uniquealleles = NULL;
	(*dat)->L_alleles = NULL;
	(*dat)->ILM = NULL;
	(*dat)->bs_ILM = NULL;
	(*dat)->pops = NULL;
	(*dat)->i_p = NULL;
	(*dat)->max_M = 0;
	return NO_ERROR;
} /* make_data */


/**
 * Free data object.
 *
 * @param dat options object
 * @return void
 */
void free_data(data *dat)
{
	int i;
	if (dat->IL)
		FREE_2ARRAY(dat->IL);
	if (dat->uniquealleles)
		FREE_1ARRAY(dat->uniquealleles);
	if (dat->L_alleles)
		FREE_2ARRAY(dat->L_alleles);
	if (dat->ILM)
		FREE_3ARRAY(dat->ILM);
	if (dat->bs_ILM)
		FREE_3ARRAY(dat->bs_ILM);
	if (dat->pops) {
		for (i=0; i<dat->numpops; i++)
			free(dat->pops[i]);
		free(dat->pops);
		dat->pops = NULL;
	}
	if (dat->i_p)
		FREE_1ARRAY(dat->i_p);

	if (dat)
		free(dat);
	dat = NULL;
} /* free_data */


/**
 * Allocate model object and initialize.
 *
 * @param mod model data object pointer reference
 * @return error status
 */
int make_model(model **mod)
{
	*mod = malloc(sizeof **mod);
	if (*mod == NULL)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, "model object");
	
	(*mod)->K = 1;
	(*mod)->pKLM = NULL;
	(*mod)->mle_pKLM = NULL;
	(*mod)->etak = NULL;
	(*mod)->mle_etak = NULL;
	(*mod)->etaik = NULL;
	(*mod)->mle_etaik = NULL;
	(*mod)->vik = NULL;
	(*mod)->diklm = NULL;
	(*mod)->I_K = NULL;
	(*mod)->IL_K = NULL;
	(*mod)->U = NULL;
	(*mod)->V = NULL;
	(*mod)->init_etaik = NULL;
	(*mod)->init_pKLM=NULL;
	(*mod)->iter1_etaik = NULL;
	(*mod)->iter2_etaik = NULL;
	(*mod)->iter1_pKLM = NULL;
	(*mod)->iter2_pKLM = NULL;
	(*mod)->current_i = 0;
	(*mod)->current_l = 0;
	(*mod)->current_k = 0;
	(*mod)->current_n = 0;
	(*mod)->etaupdate = 0;
	(*mod)->converged = 0;
	(*mod)->max_logL = -Inf;
	(*mod)->init_etak = NULL;
	(*mod)->iter1_etak = NULL;
	(*mod)->iter2_etak = NULL;
	(*mod)->em_done = 0;
	(*mod)->n_iter = 0;
	(*mod)->n_total_iter = 0;
	(*mod)->accel_abort = 0;
	return NO_ERROR;
} /* make_model */


/**
 * Allocate model.  Allocation failure leads to immediate return with error
 * code.
 *
 * [TODO] Nicer error handling.
 *
 * @param opt options object
 * @param mod model object
 * @param dat data object
 * @return error status
 */
int allocate_model(options *opt, model *mod, data *dat)
{
	MAKE_3JAGGED_ARRAY(mod->pKLM, mod->K, dat->L, dat->uniquealleles);
	if (opt->n_bootstrap && !mod->mle_pKLM)
		MAKE_3JAGGED_ARRAY(mod->mle_pKLM, mod->K, dat->L, dat->uniquealleles);

	MAKE_1ARRAY(mod->count_K, mod->K);
	
	if (opt->admixture) {
		MAKE_4JAGGED_ARRAY(mod->diklm, dat->I, mod->K, dat->L,
			dat->uniquealleles);
		MAKE_2ARRAY(mod->IL_K, dat->I*dat->ploidy, dat->L);
		MAKE_1ARRAY(mod->I_K, dat->I);
		if (opt->eta_constrained) {
			MAKE_1ARRAY(mod->etak, mod->K);
			if (opt->n_bootstrap && !mod->mle_etak)
				MAKE_1ARRAY(mod->mle_etak, mod->K);
		} else {
			MAKE_2ARRAY(mod->etaik, dat->I, mod->K);
			if (opt->n_bootstrap && !mod->mle_etaik)
				MAKE_2ARRAY(mod->mle_etaik, dat->I, mod->K);
		}
				
	} else {
		MAKE_2ARRAY(mod->vik, dat->I, mod->K);
		MAKE_1ARRAY(mod->I_K, dat->I);
		MAKE_1ARRAY(mod->etak, mod->K);
		if (opt->n_bootstrap && !mod->mle_etak)
			MAKE_1ARRAY(mod->mle_etak, mod->K);
	}
	if (opt->accel_scheme) {
		MAKE_3JAGGED_ARRAY(mod->init_pKLM, mod->K, dat->L,
			dat->uniquealleles);
		MAKE_3JAGGED_ARRAY(mod->iter1_pKLM, mod->K, dat->L,
			dat->uniquealleles);
		MAKE_3JAGGED_ARRAY(mod->iter2_pKLM, mod->K, dat->L,
			dat->uniquealleles);
		if (opt->admixture && !opt->eta_constrained) {
			MAKE_2ARRAY(mod->init_etaik, dat->I, mod->K);
			MAKE_2ARRAY(mod->iter1_etaik, dat->I, mod->K);
			MAKE_2ARRAY(mod->iter2_etaik, dat->I, mod->K);
		} else {
			MAKE_1ARRAY(mod->init_etak, mod->K);
			MAKE_1ARRAY(mod->iter1_etak, mod->K);
			MAKE_1ARRAY(mod->iter2_etak, mod->K);
		}
		MAKE_1ARRAY(mod->U, dat->max_M);
		MAKE_1ARRAY(mod->V, dat->max_M);
	}

	return NO_ERROR;
} /* allocate_model */


/**
 * Free model object.  
 *
 * @param mod model object
 * @return void
 */
void free_model_data(model *mod, options *opt)
{
	FREE_3ARRAY(mod->pKLM);
	//FREE_1ARRAY(mod->etak);
	//FREE_2ARRAY(mod->etaik);
	//FREE_2ARRAY(mod->vik);
	
	if (opt->admixture) {
		FREE_3ARRAY(mod->diklm);
		FREE_2ARRAY(mod->IL_K);
	}
	if (!opt->admixture) {
		FREE_2ARRAY(mod->vik);
		FREE_1ARRAY(mod->I_K);
	}
	//FREE_2ARRAY(mod->IL_K);
	FREE_1ARRAY(mod->count_K);
	FREE_3ARRAY(mod->init_pKLM);
	FREE_3ARRAY(mod->iter1_pKLM);
	FREE_3ARRAY(mod->iter2_pKLM);
	if (opt->admixture && !opt->eta_constrained) {
		FREE_2ARRAY(mod->init_etaik);
		FREE_2ARRAY(mod->iter1_etaik);
		FREE_2ARRAY(mod->iter2_etaik);
		FREE_2ARRAY(mod->etaik);
	} else {
		FREE_1ARRAY(mod->init_etak);
		FREE_1ARRAY(mod->iter1_etak);
		FREE_1ARRAY(mod->iter2_etak);
		FREE_1ARRAY(mod->etak);
	}
	FREE_1ARRAY(mod->U);
	FREE_1ARRAY(mod->V);
} /* free_model_data */

/**
 * Free model MLEs. MLEs are collected over multiple initializations, so the
 * memory to store them are not released until the user requests it.
 *
 * @param mod model object
 */
void free_model_mles(model *mod)
{
	if (mod->mle_pKLM)
		FREE_3ARRAY(mod->mle_pKLM);
	if (mod->mle_etak)
		FREE_1ARRAY(mod->mle_etak);
	if (mod->mle_etaik)
		FREE_2ARRAY(mod->mle_etaik);
} /* free_model_mles */

/**
 * Free model object.
 *
 * @param mod model object
 * @return void
 */
void free_model(model *mod, options *opt)
{
	if (mod) {
		free_model_data(mod,opt);
		free(mod);
	}
	mod = NULL;
} /* free_model */

/**
 * Parse command line.
 *
 * The return value contains strings that should not be freed until program
 * exit.  See print_usage() for command-line options and usage.
 *
 * @param opt options object
 * @param dat data object
 * @param argc number of command-line arguments
 * @param argv command-line arguments
 * @return error status
 */
int parse_options(options *opt, data *dat, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;

	for (i=1; i<argc; i++) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);

		/* skip to argument name */
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];

		switch (a) {
			case 'a':
				opt->admixture = 1;
				break;
			case 'b':
				opt->n_bootstrap = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->n_bootstrap < 0 || errno
					|| opt->block_relax)
					goto CMDLINE_ERROR;
				break;
			case 'c':
				opt->eta_constrained = 1;
				break;
			case 'd':
				opt->path = argv[++i];
				break;
			case 'e':
				opt->rel_error = read_double(argc, argv, ++i,
					(void *)opt);
				if (opt->rel_error < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'E':
				opt->abs_error = read_double(argc, argv, ++i,
					(void *)opt);
				if (opt->abs_error < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'f':
				opt->filename = argv[++i];
				break;
			case 'g':
				opt->adjust_step = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->adjust_step < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'h':
				fprint_usage(stdout, argv[0], (void *)opt);
				return CUSTOM_ERROR;
			case 'i':
				opt->n_init_iter = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->n_init_iter < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'I':
				opt->alleles_are_indices = 1;
				break;
			case '1':
				opt->min_K = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->min_K < 1 || errno)
					goto CMDLINE_ERROR;
			case '2':
				opt->max_K = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->max_K < 1 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'k':
				opt->max_K = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->max_K < 1 || errno)
					goto CMDLINE_ERROR;
				opt->min_K = opt->max_K;
				break;
			case 'm':
				opt->n_rand_em_init = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->n_rand_em_init == 0)
					opt->initialization_procedure = NOTHING;
				if (opt->n_rand_em_init < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'n':
				opt->n_init = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->n_init < 1 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'p':
				dat->ploidy = read_int(argc, argv, ++i,
					(void *)opt);
				if (dat->ploidy < 1 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'R':
				opt->R_format = 1;
				break;
			case 'r':
				opt->seed = read_uint(argc, argv, ++i,
					(void *)opt);
				srand(opt->seed);
				break;
			case 'x':
				opt->block_relax = 1;
				if (opt->n_bootstrap > 0 || errno)
					goto CMDLINE_ERROR;
				//if ( || errno)
				//	goto CMDLINE_ERROR;
				break;
			case 's':
				opt->accel_scheme = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->accel_scheme < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 't':
				opt->n_seconds = 60*read_uint(argc, argv, ++i,
					(void *)opt);
				break;
			case 'T':
				opt->max_iter = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->max_iter < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'v':
				opt->verbosity = VERBOSE;
				break;
			case 'P':
				opt->process = read_int(argc,argv,++i,(void*)opt);
				break;
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}

	if (opt->filename == NULL)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INVALID_CMDLINE, "You must specify the data file "
			"with command line option '-f'.  Try '-h' for help.\n");

	return NO_ERROR;

CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

/**
 * Print command-line usage.
 *
 * @param fp file stream to print on
 * @param invocation_name name of command
 * @return void
 */
void fprint_usage(FILE *fp, const char *invocation_name, void *obj)
{
	options *opt = obj ? (options *) obj : NULL;
	/* strip the "./" from the beginning of the program invocation */
	const char *prog_name = strlen(invocation_name) > 2
		&& invocation_name[0] == '.' ? &(invocation_name[2])
		: invocation_name;
	fprintf(fp, "\nNAME\n");
	fprintf(fp, "\t%s - Maximum likelihood clustering of discrete data\n",
		prog_name);
	fprintf(fp, "\nSYNOPSIS\n");
	fprintf(fp,
	"\t%s [-k <n> | -1 <n> -2 <n>] [-a -b <n> -c -d <s> -e <d> -f <d> -g <d> -h"
	"\n\t\t-i <n> -I -m <n> -n <n> -p <n> -R -s <n> -t <n> -T <d> -v -x] -f <s>\n"
	"\n\t\twhere <n> stands for integer, <s> for string",
		prog_name);
	fprintf(fp, "\nDESCRIPTION\n");
	fprintf(fp, 
	"\t%s clusters multivariate discrete data observed on a sample of\n"
	"\tindividuals using the EM algorithm.  It handles data missing at\n"
	"\trandom.  It assumes coordinates within an individual are independent.\n"
	"\tIt allows the admixture model, where each coordinate is independently\n"
	"\tdrawn from a cluster, or the mixture model, where each individual is\n"
	"\tdrawn from a cluster.\n", prog_name);
	fprintf(fp, "\nOPTIONS\n"
		/* ---------------------------------------------------------- */
		"\t-a\t"
		"Choose admixture model (default: %s).\n"
		"\t-b\t"
		"Bootstrap test of H0: K=<k>-1 vs. Ha: K=<k>, where <k> is\n"
		"\t\tgiven by -k option.  Specify number of bootstraps as\n"
		"\t\targument (default: %d).\n"
		"\t-c\t"
		"Constrain mixing proportions identical across individuals\n"
		"\t\t(only enforced with -a; default: %s).\n"
		"\t-C\t"
		"The maximum number of iterations to fit (default: %d).\n"
		"\t-d\t"
		"Directory where output files are written.\n"
		"\t-e\t"
		"Allowable log likelihood relative error for convergence\n"
		"\t\t(default: %.1e).\n"
		"\t-E\t"
		"Allowable log likelihood absolute error for convergence\n"
		"\t\t(default: %.1e).\n"
		"\t-f\t"
		"Name of data file (STRUCTURE format).\n"
		"\t-g\t"
		"Adjust step size at most this many times (default: %d)\n"
		"\t-h\t"
		"This help.\n"
		"\t-i\t"
		"Initial iterations prior to acceleration (default: %d)\n"
		"\t-I\t"
		"Alleles are indices (no sorting, etc.) (default: %s)\n"
		"\t-k\t"
		"The number of clusters to fit (default: %d).\n"
		"\t-1\t"
		"The minimum number of clusters to fit (default: %d).\n"
		"\t-2\t"
		"The maximum number of clusters to fit (default: %d).\n"
		"\t-m\t"
		"The number of Rand EM initializations, 0 to avoid Rand EM\n"
		"\t\t(default: %d).\n"
		"\t-n\t"
		"Number of initializations to run EM to convergence\n"
		"\t\t(default: %d).\n"
		"\t-p\t"
		"The ploidy (default: 2).\n"
		"\t-r\t"
		"Random number (default: %u).\n"
		"\t-R\t"
		"Data file in R format (default: %s).\n"
		"\t-s\t"
		"The acceleration scheme (default: %s).\n"
		"\t-x\t"
		"Use block relaxation algorithm (default: %s).\n"
		"\t-t\t"
		"The amount of time (in minutes) to run each k (default: %d).\n"
		"\t-v\t"
		"The verbosity level.\n",
		opt->admixture?"yes":"no",
		opt->n_bootstrap, 
		opt->eta_constrained?"yes":"no",
		opt->max_iter, opt->rel_error, opt->abs_error,
		opt->adjust_step,
		opt->n_init_iter, opt->alleles_are_indices?"yes":"no",
		opt->max_K, opt->min_K, opt->max_K, opt->n_rand_em_init,
		opt->n_init, opt->seed, opt->R_format?"yes":"no",
		accel_method_names[opt->accel_scheme],
		opt->block_relax?"yes":"no", opt->n_seconds
	);
} /* fprintf_usage */
