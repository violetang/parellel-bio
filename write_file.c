/**
 * @file write_file.c
 * @author Arun Sethuraman
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Wed Dec  5 20:44:34 CST 2012
 *
 * Functions to write results to output files.
 */

#include "multiclust.h"

/**
 * Write data to file.  Presumably, you call this function to write
 * simulated data, and simulated data is in _data::ILM, so set last
 * argument non-zero.
 *
 * @param opt options ojbect
 * @param dat data object
 * @param use_ILM write from dat->ILM (not dat->IL)
 * @return error status
 */
int write_data(options *opt, data *dat, int use_ILM)
{
	int i, j, l, m, msum, r, r2;
	char *outfile = NULL;
	int len = strlen(opt->path) + 7;
	FILE *fp;

	r2 = r = rand();
	while (r2) {
		len++;
		r2 /= 10;
	}

	MAKE_1ARRAY(outfile, len);

	sprintf(outfile, "%sbs%d.str", opt->path, r);

	if ((fp = fopen(outfile, "w")) == NULL) {
		FREE_1ARRAY(outfile);
		message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_OPEN_ERROR, outfile);
	}

	/* header line */
	for (l=0; l<dat->L; l++)
		fprintf(fp, "%sloc%d", l?" ":"", l+1);
	fprintf(fp, "\n");

	if (!use_ILM) {
		for (i=0; i<dat->I; i++)
			for (j=1; j<=dat->ploidy; j++) {
				for (l=0; l<dat->L; l++)
					fprintf(fp, " %d", dat->IL[i][l]);
				fprintf(fp, "\n");
			}
	} else {

		for (i=0; i<dat->I; i++) {			/* individual */
			for (j=1; j<=dat->ploidy; j++) {	/* haplotype */
				fprintf(fp, "%s %s", dat->idv[i].name,
					dat->pops[dat->idv[i].locale]);
				for (l=0; l<dat->L; l++) {	/* locus */
					/* skip ahead to next allele to print */
					m = 0;
					msum = dat->ILM[i][l][m];
					while (msum < j)
						msum += dat->ILM[i][l][++m];
					if (dat->L_alleles)
						fprintf(fp, " %d", dat->L_alleles[l][m]);
					else
						fprintf(fp, " %d", m);
				}
				fprintf(fp, "\n");
			}
		}
	}
	
	fclose(fp);

	FREE_1ARRAY(outfile);

	return NO_ERROR;
} /* End of write_data(). */


/**
 * Write maximum likelihood parameter estimates to file.
 *
 * @param outfile name of output file
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int write_file(char *outfile, options *opt, data *dat, model *mod)
{
	FILE *fp;
	int i, k, l, m;

	if ((fp = fopen(outfile, "w")) == NULL)
		return message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, FILE_OPEN_ERROR, outfile);

	/* mixing proportions eta */
	if (!opt->admixture || opt->eta_constrained) {
		/* etak */
		for(k = 0; k < mod->K; k++)
			fprintf(fp, "%f ", mod->etak[k]);
		fprintf(fp, "\n\n");
	} else {
		/* etaik */
		for (i = 0; i < dat->I; i++) {
			for (k = 0; k < mod->K; k++)
				fprintf(fp, "%f ", mod->etaik[i][k]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n");
	}

	/* allele probabilities pKLM */
	for (k = 0; k < mod->K; k++) {
		for (l = 0; l < dat->L; l++) {
			for (m = 0; m < dat->uniquealleles[l]; m++)
				fprintf(fp, "%f ", mod->pKLM[k][l][m]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	return NO_ERROR;
} /* End of write_file(). */

/**
 * Write detailed output about fitted admixture model.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int write_file_detail(options *opt, data *dat, model *mod)
{
	FILE *fp;
	int i, k, l, m, m_start;
	int len;		/* length of outfile name */
	double p;		/* number of parameters */
	char *outfile = NULL;	/* outfile name */

	/* count number of parameters */
	if (!opt->admixture || opt->eta_constrained) 
		p = (mod->K - 1);
	else
		p = dat->I * (mod->K - 1);
	for (l = 0; l < dat->L; l++) {
		/* was BUG: forgot to not count missing allele */
		m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING
			? 1 : 0;
		p += (dat->uniquealleles[l] - 1 - m_start) * mod->K;
	}

	/* compute length of output filename */
	len = strlen(opt->path) + strlen(opt->filename) + 14;
	/* add number of digits for mod->K */
	k = mod->K;
	while (k) {
		len++;
		k /= 10;
	}

	MAKE_1ARRAY(outfile, len);

	/* print maximum log likelihood, AIC, BIC, and partition counts */
	sprintf(outfile, "%s%s_out_K=%d.txt", opt->path, opt->filename, mod->K);
	if ((fp = fopen(outfile, "w")) == NULL) {
		FREE_1ARRAY(outfile);
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_OPEN_ERROR, outfile);
	}

	fprintf(fp, "logL = %f\n", mod->logL);
	fprintf(fp, "AIC = %f\n", -2 * mod->logL + 2 * p);
	fprintf(fp, "BIC = %f\n\n", -2 * mod->logL + p * log((double) dat->I));
	fprintf(fp, "count.K\n");
	for (k = 0; k < mod->K; k++)
		fprintf(fp, "%d ", mod->count_K[k]);
	fprintf(fp, "\n\n");
	fclose(fp);

	/* print etaik */
	if (!opt->admixture || opt->eta_constrained) {
		sprintf(outfile, "%s%s_etak_K=%d.txt", opt->path, opt->filename,
			mod->K);
		if ((fp = fopen(outfile, "w")) == NULL)
			return message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, FILE_OPEN_ERROR, outfile);

		fprintf(fp, "i\tk\tetak\n");
		for (k = 0; k < mod->K; k++)
			fprintf(fp, "%d\t%f\n", k, mod->etak[k]);
		fprintf(fp, "\n");
	} else {
		sprintf(outfile, "%s%s_etaik_K=%d.txt", opt->path, opt->filename,
			mod->K);
		if ((fp = fopen(outfile, "w")) == NULL)
			return message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, FILE_OPEN_ERROR, outfile);

		fprintf(fp, "i\tk\tetaik\n");
		for (i = 0; i < dat->I; i++)
			for (k = 0; k < mod->K; k++)
				fprintf(fp, "%d\t%d\t%f\n", i, k, mod->etaik[i][k]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	/* print KLM */
	sprintf(outfile, "%s%s_klm_K=%d.txt", opt->path, opt->filename, mod->K);
	if ((fp = fopen(outfile, "w")) == NULL)
		return message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, FILE_OPEN_ERROR, outfile);

	fprintf(fp, "k\tl\tm\tKLM\n");
	for (k = 0; k < mod->K; k++)
		for (l = 0; l < dat->L; l++)
			for (m = 0; m < dat->uniquealleles[l]; m++)
				fprintf(fp, "%d\t%d\t%d\t%f\n", k, l, m,
					mod->pKLM[k][l][m]);
	fprintf(fp, "\n");
	fclose(fp);

	FREE_1ARRAY(outfile);

	return NO_ERROR;
} /* End of write_file_admixture_detail(). */


/**
 * Assign each individual to a posteriori most likely subpopulation.  The
 * admixture model allows alleles to arise from different subpopulations, so
 * the assignment here is based on which subpopulation is the a posterior most
 * probable source of the most alleles.
 *
 * was BUG: was doing a priori assignments based on estimated etaik
 *
 * @param dat data object
 * @param mod model object
 * @return void
 */
void partition_admixture(data *dat, model *mod)
{
	int i, k, l, m, I_K;
	double map;
	double dik[mod->K];

	for (k = 0; k < mod->K; k++)
		mod->count_K[k] = 0;

	for (i = 0; i < dat->I; i++) {
		/* compute expected number of alleles in subpopulation k */
		for (k = 0; k < mod->K; k++)
			dik[k] = 0;
		for (l = 0; l < dat->L; l++)
			for (m = 0; m < dat->uniquealleles[l]; m++)
				for (k = 0; k < mod->K; k++)
					dik[k] += mod->diklm[i][k][l][m];

		/* assign individual to subpopulation source of most alleles */
		I_K = 0;
		map = dik[0];
		for (k = 1; k < mod->K; k++)
			if (dik[k] > map) {
				map = dik[k];
				I_K = k;
			}

		/* assign to subpopulation I_K; increment I_K subpop. cnt */
		mod->I_K[i] = I_K;
		mod->count_K[I_K]++;
	}
} /* End of partition_admixture(). */


/**
 * Print POPQ file.  This function prints out the POPQ file, as used by
 * DISTRUCT/CLUMPP, under the admixture model.
 *
 * [TODO] Sorry, error handling is a really awkward aspect of array.h.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
#undef MAKE_1ARRAY
#define MAKE_1ARRAY MAKE_1ARRAY_ERR
int popq_admix(options *opt, data *dat, model *mod)
{
	int i, k, l, m, n;
	int len;
	char *outfile = NULL;
	FILE *fp;
	double **vik_p = NULL;
	int err = NO_ERROR;

	len = strlen(opt->path) + strlen(opt->filename) + 18;
	k = mod->K;
	while (k) {
		len++;
		k /= 10;
	}

	MAKE_1ARRAY(outfile, len);
	if ((err = errno))
		goto FREE_AND_RETURN;

	/* open file to write */
	sprintf(outfile, "%s%s_admix_popq_%d.popq", opt->path, opt->filename,
		mod->K);
	if ((fp = fopen(outfile, "w")) == NULL) {
		err = message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_OPEN_ERROR, outfile);
		goto FREE_AND_RETURN;
	}

	/* compute post. prob. that locale n allele is in pop. k */
	MAKE_2ARRAY(vik_p, dat->numpops, mod->K);
	if ((err = errno))
		goto FREE_AND_RETURN;
	for (n = 0; n < dat->numpops; n++)
		for (k = 0; k < mod->K; k++)
			vik_p[n][k] = 0;

	for (k = 0; k < mod->K; k++) {
		for (i = 0; i < dat->I; i++)
			for (l = 0; l < dat->L; l++) {
				/* was BUG: skipped missing data, but
				 * - missing data has info, and
				 * - missing alleles used to normalize below
				 */
				for (m=0; m<dat->uniquealleles[l]; m++)
					vik_p[dat->idv[i].locale][k]
						+= mod->diklm[i][k][l][m];
			}
		for (n = 0; n < dat->numpops; n++)
			vik_p[n][k] = vik_p[n][k]
				/(dat->ploidy*dat->L*dat->i_p[n]);
	}

	for (n = 0; n < dat->numpops; n++) {
		fprintf(fp, "%s:\t", dat->pops[n]);
		for (k = 0; k < mod->K; k++)
			fprintf(fp, "%lf\t", vik_p[n][k]);
		fprintf(fp, "%d\n", dat->i_p[n]);
	}

	fclose(fp);

FREE_AND_RETURN:
	FREE_2ARRAY(vik_p);
	FREE_1ARRAY(outfile);

	return err;
} /* end of popq_admix(). */
#undef MAKE_1ARRAY
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN


/**
 * Print out INDIVQ file, as used by DISTRUCT/CLUMPP, under the admixture model.
 *
 * [TODO] Sorry, error handling is a really awkward aspect of array.h.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
#undef MAKE_1ARRAY
#define MAKE_1ARRAY MAKE_1ARRAY_ERR
int indivq_admix(options *opt, data *dat, model *mod)
{

	int i, k, l, m;
	int len;
	char *outfile = NULL;
	double **vik = NULL;
	FILE *fp;
	int err = NO_ERROR;

	len = strlen(opt->path) + strlen(opt->filename) + 22;
	k = mod->K;
	while (k) {
		len++;
		k /= 10;
	}

	MAKE_1ARRAY(outfile, len);
	if ((err = errno))
		goto FREE_AND_RETURN;

	sprintf(outfile, "%s%s_admix_indivq_%d.indivq", opt->path, opt->filename, mod->K);
	if ((fp = fopen(outfile, "w")) == NULL) {
		err = message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_OPEN_ERROR, outfile);
		goto FREE_AND_RETURN;
	}

	/* compute probability random allele from i from subpopn k */
	MAKE_2ARRAY(vik, dat->I, mod->K);
	if ((err = errno))
		goto FREE_AND_RETURN;

	for (k = 0; k < mod->K; k++)
		for (i = 0; i < dat->I; i++)
			vik[i][k] = 0;

	for (k = 0; k < mod->K; k++)
		for (i = 0; i < dat->I; i++) {
			for (l = 0; l < dat->L; l++)
				/* was BUG ignored missing as in popq_admix() */
				for (m = 0; m < dat->uniquealleles[l]; m++)
					vik[i][k] += mod->diklm[i][k][l][m];
			vik[i][k] = vik[i][k]/(dat->ploidy*dat->L);
		}

	for (i = 0; i < dat->I; i++) {

		fprintf(fp, "%d\t%s\t(x)\t%s\t:", i, dat->idv[i].name,
			dat->pops[dat->idv[i].locale]);
		for (k = 0; k < mod->K; k++)
			fprintf(fp, "\t%f", vik[i][k]);
		fprintf(fp, "\n");
	}

	fclose(fp);

FREE_AND_RETURN:
	FREE_2ARRAY(vik);
	FREE_1ARRAY(outfile);

	return err;
} /* end of indivq_admix(). */
#undef MAKE_1ARRAY
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN


/**
 * Assign individuals to subpopulations.  Assign individuals to subpopulation
 * based on Maximum A Posteriori (MAP) estimate, assuming flat prior.
 *
 * @param dat data object
 * @param mod model object
 * @return void
 */
void partition_mixture(data *dat, model *mod)
{
	int i, k;
	double tmp_vik;

	for (k = 0; k < mod->K; k++)
		mod->count_K[k] = 0;

	for (i = 0; i < dat->I; i++) {
		mod->I_K[i] = 0;
		tmp_vik = mod->vik[i][0];
		for (k = 1; k < mod->K; k++)
			if (mod->vik[i][k] > tmp_vik) {
				tmp_vik = mod->vik[i][k];
				mod->I_K[i] = k;
			}
		mod->count_K[mod->I_K[i]]++;
	}
} /* End of partition_mixture(). */


/**
 * Print a popq file.  Print data in a popq-style file under the mixture model
 * as required by DISTRUCT/CLUMPP with sufficient statistics.
 *
 * [TODO] Sorry, error handling is a really awkward aspect of array.h.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
#undef MAKE_1ARRAY
#define MAKE_1ARRAY MAKE_1ARRAY_ERR
int popq_mix(options *opt, data *dat, model *mod)
{
	int i, k, n, len;
	char *outfile = NULL;
	double **vik_p = NULL;
	FILE *fp;
	int err = NO_ERROR;

	/* compute length of output filename */
	len = strlen(opt->path) + strlen(opt->filename) + 14;

	MAKE_1ARRAY(outfile, len);
	if ((err = errno))
		goto FREE_AND_RETURN;

	/* write popq file */
	sprintf(outfile, "%s%smix_popq.popq", opt->path, opt->filename);
	if ((fp = fopen(outfile, "w")) == NULL) {
		err = message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, FILE_OPEN_ERROR, outfile);
		goto FREE_AND_RETURN;
	}

	MAKE_2ARRAY(vik_p, dat->numpops, mod->K);
	if ((err = errno))
		goto FREE_AND_RETURN;

	for (n=0; n<dat->numpops; n++)
		for (k = 0; k < mod->K; k++)
			vik_p[n][k] = 0;

	for (k = 0; k < mod->K; k++)		/* population */
		/* was BUG: access past end of array */
		for (i=0; i<dat->I; i++)	/* individual */
			vik_p[dat->idv[i].locale][k] += mod->vik[i][k];

	/* normalize by the number of individuals in each locale */
	for (n = 0; n < dat->numpops; n++)
		for (k = 0; k < mod->K; k++)
			vik_p[n][k] = vik_p[n][k]/dat->i_p[n];

	for (n = 0; n < dat->numpops; n++) {
		fprintf(fp, "%s:\t", dat->pops[n]);
		for (k = 0; k < mod->K; k++)
			fprintf(fp, "%lf\t", vik_p[n][k]);
		fprintf(fp, "%d\n", dat->i_p[n]);
	}

	fclose(fp);

FREE_AND_RETURN:
	FREE_2ARRAY(vik_p);
	FREE_1ARRAY(outfile);

	return err;
} /* end of popq_mix() */
#undef MAKE_1ARRAY
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN


/**
 * Print results in INDIVQ suitable for DISTRUCT/CLUMPP.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 * @return error status
 */
int indivq_mix(options *opt, data *dat, model *mod)
{
	int i, k, len;
	char *outfile = NULL;
	FILE *fp;

	len = strlen(opt->path) + strlen(opt->filename) + 18;

	MAKE_1ARRAY(outfile, len);

	sprintf(outfile, "%s%smix_indivq.indivq", opt->path, opt->filename);
	if ((fp = fopen(outfile, "w")) == NULL) {
		FREE_1ARRAY(outfile);
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_OPEN_ERROR, outfile);
	}

	for (i = 0; i < dat->I; i++) {

		fprintf(fp, "%d\t%s\t(x)\t%s\t:", i, dat->idv[i].name,
			dat->pops[dat->idv[i].locale]);
		for (k = 0; k < mod->K; k++)
			fprintf(fp, "\t%f", mod->vik[i][k]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	FREE_1ARRAY(outfile);

	return NO_ERROR;
} /* end of indivq_mix() */
