#' Regression models for local mutations: Evaluating differential mutation rates across classes of sites
#'
#' RM2() uses negative binomial regression to evaluate local mutation frequencies and processes between stacked sites of the same class to flanking control regions
#' @param maf Data frame of mutations (prepared by get_mut_trinuc_strand) containiing the following information:
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the mutation in base 1 coordinates}
#'     \item{end}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#'     \item{mut_trinuc}{trinucleotide context - where middle is C or T - with alternate allele}
#'     \item{mut_strand}{character indicating Watson (w) or Crick (c)}
#'     \item{ref_alt}{character indicating single-base substitution}
#' }
#' @param sites Data frame of site coordinates
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the mutation in base 0 coordinates}
#'     \item{end}{the end position of the mutation in base 0 coordinates}
#' }
#' @param mut_class_column Character corresponding to column of mutation classes for grouped analysis
#' @param cofactor_column Character corresponding to column of binary cofactors
#' @param window_size Integer indicating the half-width of sites and flanking regions (added to left and right for full width). (default 100)
#' @param n_min_mut Integer indicating the minimum number of mutations required to perform analysis (default 100)
#' @param n_bin Integer indicating the number of megabase bins to use (default 10)
#'
#' @return Data frame containing the regression estimates and likelihood ratio test output with the following columns: mut_type,
#' pp, this_coef, obs_mut, exp_mut, exp_mut_hi, exp_mut_lo, fc, n_sites_tested
#' \describe{
#'     \item{mut_type}{A string identifying the mutation class}
#'     \item{pp}{The p-value from the likelihood ratio test}
#'     \item{this_coef}{The coefficient from is_site}
#'     \item{obs_mut}{The total number of observed mutations of that class}
#'     \item{exp_mut}{The expected number of mutations determined by the model}
#'     \item{exp_mut_lo}{Lower bound of 95% confidence interval}
#'     \item{exp_mut_hi}{Upper bound of 95% confidence interval}
#'     \item{fc}{Observed mutations divided by expected mutations}
#'     \item{pp_cofac}{The p-value from the likelihood ratio test of site:cofactor interaction}
#'     \item{this_coef_cofac}{The coefficient from the site:cofactor interaction term}
#'     \item{n_sites_tested}{The number of sites that were tested - all sites if no downsampling}
#' }
#' @export
RM2 = function(maf, sites, mut_class_columns = NA, cofactor_column = NA,
			window_size = 100, n_min_mut = 100, n_bin = 10) {

	# return this boilerplate if errors caught
	empty_res = data.frame(mut_type = NA, pp = NA, this_coef = NA, obs_mut = NA, 
          			exp_mut = NA, exp_mut_lo = NA, exp_mut_hi = NA, fc = NA, 
          			pp_cofac = NA, this_coef_cofac = NA, n_sites_tested = NA,
          			stringsAsFactors = FALSE)
	
	prepared_sites = NULL
  if (is.na(cofactor_column)) {
  		prepared_sites = try(.prepare_bins_of_sites(sites, window_size, maf, 
  						                n_bin = n_bin, n_min_mut = n_min_mut), silent = TRUE)
  		if (any(class(prepared_sites) == "try-error")) {
  			
    			# check if error is one of permitted errors
    			decide = .shall_it_pass(prepared_sites)
    			if (decide) {
    			    return(empty_res)
    			} else {
    				  stop(prepared_sites)
    			}
  		}
  	} else {
		split_maf = split(maf, maf[,cofactor_column])
		prepared_sites = try(lapply(names(split_maf), function(cof) 
				                    .prepare_bins_of_sites(sites, window_size, split_maf[[cof]], 
					                	                        n_bin = n_bin, n_min_mut = n_min_mut)), silent = TRUE)
		if (any(class(prepared_sites) == "try-error")) {
  			decide = .shall_it_pass(prepared_sites)
  			if (decide) {
  				  return(empty_res)
  			} else {
  				  stop(prepared_sites)
  			}
		}
		names(prepared_sites) = try(names(split_maf), silent = TRUE)
		rm(split_maf)
	}

	maf_list = do.call(c, lapply(mut_class_columns, function(mc) {
  		if (is.na(mc)) {
  			  list("total_muts__total_muts" = maf)
  		} else {
    			res = split(maf, maf[,mc])
    			names(res) = paste(mc, names(res), sep="__")
    			res
  		}
	}))

	# maf1_title = "total_muts__total_muts"
	res = do.call(rbind, lapply(names(maf_list), function(mc) {
			maf1 = maf_list[[mc]]

  		dfr = NULL
  		if (is.na(cofactor_column)) {
  			  dfr = try(.maf_to_dfr2(maf1, prepared_sites), silent = TRUE)
  		} else {
    			split_maf1 = split(maf1, maf1[,cofactor_column])
    			dfr = try(do.call(rbind, lapply(names(split_maf1), function(cof) {
            					split_dfr = .maf_to_dfr2(split_maf1[[cof]], prepared_sites[[cof]])
            					split_dfr = data.frame(split_dfr, cof, stringsAsFactors = FALSE)
            					split_dfr
    			})), silent = TRUE)
  		}

		if (any(class(dfr) == "try-error")) {
  		  cat(dfr, "\n")
  			dfr = NULL
		}
    rm(maf1)

		this_res = try(.test_NB_model(dfr, mc, n_min_mut = n_min_mut,
						        test_cofactor = !is.na(cofactor_column)), silent = TRUE)

				if (any(class(this_res) == "try-error")) {
  				# check if error is one of permitted errors
    			decide = .shall_it_pass(this_res)
    			if (decide) {
      				empty_res$mut_type = mc
      				this_res = empty_res
    			} else {
      				stop(this_res)
    			}
				}
				this_res$n_sites_tested = nrow(sites)
				this_res
			}))

			rownames(res) = NULL
			res
}


#' Determines whether error is permitted
#'
#' @param this_err Character error message
#'
#' @return Boolean indicating whether the error is permitted or not
.shall_it_pass = function (this_err) {

	this_err1 = gsub("Error : (.+)\n$", "\\1", this_err)	
	permitted_err = c(
              			"Too few mutations in sites+flanks.",
              			"glm.nb failed.",
              			"Too few mutations to produce bins.")
	this_err1 %in% permitted_err
}


#' Evaluating differential mutation rates across classes of sites with downsampling
#'
#' RM2_downsample() is a wrapper that first performs downsampling then calls RM2. The median index is selected by p-value for total mutations (mut_class_columns=NA) and the corresponding values are returned.
#' @param maf Data frame of mutations
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the mutation in base 1 coordinates}
#'     \item{end}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#'     \item{mut_trinuc}{trinucleotide context - where middle is C or T - with alternate allele}
#'     \item{mut_strand}{character indicating Watson (w) or Crick (c)}
#'     \item{ref_alt}{character indicating single-base substitution}
#' }
#' @param sites Data frame
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the mutation in base 0 coordinates}
#'     \item{end}{the end position of the mutation in base 0 coordinates}
#' }
#' @param mut_class_column Character corresponding to column of mutation classes for grouped analysis
#' @param cofactor_column Character corresponding to column of cofactors
#' @param window_size Integer indicating the half-width of sites and flanking regions (added to left and right for full width) (default 100)
#' @param n_min_mut Integer indicating the minimum number of mutations required to perform analysis (default 100)
#' @param n_bin Integer indicating the number of megabase bins to use (default 10)
#' @param n_sites_sampled Integer indicating the number of sites to sample
#' @param n_iterations Integer indicating how many times to repeat the sampling procedure (default 100)
#'
#' @return Data frame containing the regression estimates and likelihood ratio test output with the following columns: mut_type,
#' pp, this_coef, obs_mut, exp_mut, exp_mut_hi, exp_mut_lo, fc, n_sites_tested
#' \describe{
#'     \item{mut_type}{A string identifying the mutation class}
#'     \item{pp}{The p-value from the likelihood ratio test}
#'     \item{this_coef}{The coefficient from is_site}
#'     \item{obs_mut}{The total number of observed mutations of that class}
#'     \item{exp_mut}{The expected number of mutations determined by the model}
#'     \item{exp_mut_lo}{Lower bound of 95% confidence interval}
#'     \item{exp_mut_hi}{Upper bound of 95% confidence interval}
#'     \item{fc}{Observed mutations divided by expected mutations}
#'     \item{pp_cofac}{The p-value from the likelihood ratio test of site:cofactor interaction}
#'     \item{this_coef_cofac}{The coefficient from the site:cofactor interaction term}
#'     \item{n_sites_tested}{The number of sites that were tested - all sites if no downsampling}
#' }
#' @export
RM2_downsample = function(maf, sites, mut_class_columns = NA, cofactor_column = NA, 
                              window_size = 100, n_min_mut = 100, n_bin = 10, 
                              n_sites_sampled, n_iterations = 100) {
  
  if (nrow(sites) < n_sites_sampled) {
    res1 = RM2(maf, sites, mut_class_columns = mut_class_columns, 
                   cofactor_column = cofactor_column, window_size = window_size, 
                   n_min_mut = n_min_mut, n_bin = n_bin)
    return(res1)
  }
  
  sampled_res = lapply(1:n_iterations, function(i) {
    cat(".")
    sites1 = sites[sample(1:nrow(sites), n_sites_sampled),]
    res1 = RM2(maf, sites1, mut_class_columns = mut_class_columns, 
                   cofactor_column = cofactor_column, window_size = window_size, 
                   n_min_mut = n_min_mut, n_bin = n_bin)
  })
  
  # select representative by median p-value
  select_index = .which_median(sapply(sampled_res, function(x) 
    x[x$mut_type == "total_muts__total_muts","pp"]))
  
  sampled_res[[select_index]]
}



#' Runs negative binomial regression and likelihood ratio test
#'
#' @param dfr Data frame generated from maf_to_dfr that contains mutation and trinucleotide counts for each quadnucleotide context + indel
#' @param mut_type Character corresponding to total mutations or name of subclass
#' @param n_min_mut Integer (default 100) representing the minimum number of mutations required to run test
#' @param test_cofactor Boolean indicating whether to analyze cofactor
#'
#' @return Data frame containing the results of the regression with the following columns: mut_type,
#' pp, this_coef, obs_mut, exp_mut, exp_mut_lo, exp_mut_hi, fc, pp_cofac, this_coef_cofac
#' \describe{
#'     \item{mut_type}{A string identifying the mutation class}
#'     \item{pp}{The p-value from the likelihood ratio test}
#'     \item{this_coef}{The coefficient from is_site}
#'     \item{obs_mut}{The total number of observed mutations of that class}
#'     \item{exp_mut}{The expected number of mutations determined by the model}
#'     \item{exp_mut_lo}{Lower bound of 95% confidence interval}
#'     \item{exp_mut_hi}{Upper bound of 95% confidence interval}
#'     \item{fc}{Observed mutations divided by expected mutations}
#'     \item{pp_cofac}{The p-value from the likelihood ratio test of site:cofactor interaction}
#'     \item{this_coef_cofac}{The coefficient from the site:cofactor interaction term}
#' }
.test_NB_model = function(dfr, mut_type, n_min_mut = 100, test_cofactor = F) {

	if (test_cofactor & length(unique(dfr$cof)) != 2) {
      	stop("Incorrect number of cofactor values.", call. = FALSE)
	}

	if (sum(dfr$muts_count) < n_min_mut) {
    		stop("Too few mutations in sites+flanks.", call. = FALSE)
	}

	silent_quadnucs = names(which(c(by(dfr$muts_count, dfr$quadnuc, sum))==0))
	dfr = dfr[!dfr$quadnuc %in% silent_quadnucs,, drop=F ]

	# mutation rate and count of positions will be log'd for model stability
	dfr$log_mut_rate = log1p(dfr$mut_rate)
	dfr$log_posi_count = log(dfr$posi_count)

	this_formula = "muts_count ~ log_mut_rate + is_site + offset(log_posi_count)"

	# adjust formula if only one quad-nuc context of mutations is available; most commonly indel models
	if (length(unique(dfr$quadnuc)) > 1)	{
	    this_formula = paste(this_formula, "+ quadnuc")
	}

	if (test_cofactor)	{
		  this_formula = paste(this_formula, "+ cof")
	}

	h1 = try(MASS::glm.nb(stats::as.formula(this_formula), data=dfr), silent=T)

	if (any(class(h1)=="try-error")) {
			stop("glm.nb failed.", call. = FALSE)
	}

	h0 = try(stats::update(h1, . ~ . - is_site, data=dfr), silent=T)

	if (any(class(h0)=="try-error")) {
		  stop("glm.nb failed.", call. = FALSE)
	}

	# "manual" calculation is more precise, thanks to lower.tail.
	# if using 1-p(), again limited to 10-16.
	# anova(h1) provides wald test that is less conservative
	pp = stats::pchisq(-(h0$twologlik - h1$twologlik), 1, lower.tail=F)
	this_coef = stats::coef(h1)[['is_siteTRUE']]

	# compute predicted mutation counts with ranges
	h0_mut = .predict_muts_ranges_probs(h0, which(dfr$is_site))
	h1_mut = .predict_muts_ranges_probs(h1, which(dfr$is_site))
	real_obs_mut = sum(dfr[dfr$is_site,"muts_count"])
	exp_mut = h0_mut[['mid']]
	exp_mut_lo = h0_mut[['lo']]
	exp_mut_hi = h0_mut[['hi']]
	fc = h1_mut[['mid']] / h0_mut[['mid']]
	
	pp_cofac = this_coef_cofac = NA
	if (test_cofactor) {

		h2 = try(stats::update(h1, . ~ . + cof:is_site, data=dfr), silent=T)
		
		if (any(class(h2) == "try-error")) {
	    	stop("glm.nb failed.", call. = FALSE)
		}

		pp_cofac = stats::pchisq(-(h1$twologlik - h2$twologlik), 1, lower.tail=F)
		# coefficient is computed with first level of factor as reference
		this_coef_cofac_name = paste0("is_siteTRUE:cof", sort(unique(dfr$cof))[2])
		this_coef_cofac = stats::coef(h2)[[this_coef_cofac_name]]
	}

	data.frame(mut_type, pp, this_coef, obs_mut=real_obs_mut,
    		exp_mut, exp_mut_lo, exp_mut_hi, fc, pp_cofac, this_coef_cofac,
    		stringsAsFactors=F)
}



#' Extract expected mutations count estimates by probability sampling
#'
#' @param hyp List containing glm.nb regression outputs
#' @param select_positions Numeric vector containing indices of site positions
#'
#' @return Vector of mid, lo and hi mutation count predictions
.predict_muts_ranges_probs = function(hyp, select_positions) {
  hyp_predict_probs = replicate(1000, 
                          sum(MASS::rnegbin(stats::fitted(hyp), theta = (summary(hyp))$theta)[select_positions]))
	mid = stats::median(hyp_predict_probs)
	lo = stats::quantile(hyp_predict_probs, 0.025, na.rm=T)[[1]]
	hi = stats::quantile(hyp_predict_probs, 0.975, na.rm=T)[[1]]

	c(mid=mid, lo=lo, hi=hi)
}



#' Selects median value with priority given to larger value
#'
#' @param x Numeric vector of p-values
#'
#' @return Position of median value
.which_median = function(x) {
  if (all(is.na(x))) { return(1) }
    mid_pos = ceiling(length(x)/2) 
    index = which(x == sort(x)[mid_pos])[1]
  if (is.na(index)) {
    index = mid_pos
  }
  index
}


