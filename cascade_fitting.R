# Filename: cascade_fitting.R
#
# Authors:
# 	Linh Huynh - huynh@ucdavis.edu
# 	Jason Youn - jyoun@ucdavis.edu
#
# Description:
#
# To-do:
#

#################
# initial setup #
#################

# setup directories and file names
project_parent_dir <- getwd()
project_data_dir <- file.path(project_parent_dir, "data")
pTET_data_filepath <- file.path(project_data_dir, "pTET.dat")
pBAD_data_filepath <- file.path(project_data_dir, "pBAD.dat")
pBAD_pTET_data_filepath <- file.path(project_data_dir, "pBAD_pTET.dat")

message("pTET data filepath: ", pTET_data_filepath)
message("pBAD data filepath: ", pBAD_data_filepath)
message("pBAD_pTET data filepath: ", pBAD_pTET_data_filepath)
cat("\n")

# read data
pTET_mat <- read.table(pTET_data_filepath, header=TRUE)
pBAD_mat <- read.table(pBAD_data_filepath, header=TRUE)
pBAD_pTET_mat <- read.table(pBAD_pTET_data_filepath, header=TRUE)

message(sprintf("pTET size: (%d, %d)", nrow(pTET_mat), ncol(pTET_mat)))
message(sprintf("pBAD size: (%d, %d)", nrow(pBAD_mat), ncol(pBAD_mat)))
message(sprintf("pBAD_pTET size: (%d, %d)", nrow(pBAD_pTET_mat), ncol(pBAD_pTET_mat)))
cat("\n")


###########################
# fit parameters for pTET #
###########################

# define the model where each pTET mutant has 4 parameters:
#	a: basal expression
#	b: strength
#	k: affinity
#	n: cooperativity
tetR_const = 1
k_aTc = 10
n_aTc = 4
tetR_sigmoid <- y ~ a + b / (1 + (tetR_const / ((1 + (aTc / k_aTc) ^ n_aTc) * k)) ^ n)

# create empty zero matrix to populate the parameters below
pTET_parameter <- matrix(0, nrow=4, ncol=(ncol(pTET_mat) - 1)) # (4, 18)

# create empty zero matrix to populate both actual and prediction
pTET_pred <- matrix(0, nrow=nrow(pTET_mat), ncol=ncol(pTET_mat))
pTET_pred[, 1] = pTET_mat[, 1]

# fit the model and find parameters
message("Fitting the model for pTET...")
cat("\n")

for (mutant in (2:ncol(pTET_mat))) {
	message("Fitting mutant ", colnames(pTET_mat)[mutant])

	# construct tetR_data to fit
	ydat <- pTET_mat[, mutant] # mutant, y-axis
	tdat <- pTET_mat[, 1] # aTc, x-axis
	tetR_data <- data.frame(y=ydat, aTc=tdat)

	# define starting point and range of the parameters to feed into nls
	a1 = min(ydat) / 10
	a2 = min(ydat) * 2
	b1 = max(ydat) / 100
	b2 = max(ydat) * 2
	lower_bound <- c(a=a1, b=b1, k=0.01, n=1)
	upper_bound <- c(a=a2, b=b2, k=1, n=4)
	init_val <- c(a=a1, b=b1, k=0.1, n=2)
	obj_w <- 1 / ydat

	# use NLS (non-linear least square) to fit 4 parameters of each mutant with the experimental data
	anlsb1 <- try(nls(tetR_sigmoid, start=init_val, lower=lower_bound, upper=upper_bound, data=tetR_data, algorithm="port", weights=obj_w))
	tetR_data$predict <- predict(anlsb1, interval="predict") # append prediction

	message("tetR_data:")
	print(tetR_data)

	# extract model coefficients
	coefficients = coef(anlsb1)
	pTET_parameter[, mutant-1] = coefficients # populate the coefficients to pTET_parameter

	# save predictions
	pTET_pred[, mutant] = tetR_data$predict

	cat("\n")
}

# By now, all parameter values were found and stored in the matrix pTET_parameter.
# In the matrix pTET_parameter, rows represent parameters and columns represent mutants.
message("pTET_parameter:")
print(pTET_parameter)

# save the predictions
write.table(pTET_pred, file=file.path(project_data_dir, "pTET_pred.txt"), row.names=FALSE, col.names=colnames(pTET_mat), sep="\t")


###########################
# fit parameters for pBAD #
###########################

# define the model where each pDAT mutant also has same 4 parameters:
#	a: basal expression
#	b: strength
#	k: affinity
#	n: cooperativity
araC_const = 1
k_Lara = 0.00001
n_Lara = 2
AraC_sigmoid <- y ~ a + b / (1 + (araC_const / ((1 + (Lara / k_Lara) ^ n_Lara) * k)) ^ n)

# create empty zero matrix to populate the parameters below
pBAD_parameter <- matrix(0, nrow=4, ncol=ncol(pBAD_mat) - 1)

# create empty zero matrix to populate both actual and prediction
pBAD_pred <- matrix(0, nrow=nrow(pBAD_mat), ncol=ncol(pBAD_mat))
pBAD_pred[, 1] = pBAD_mat[, 1]

# fit the model and find parameters
message("Fitting the model for pBAD...")
cat("\n")

for (mutant in (2:ncol(pBAD_mat))) {
	message("Fitting mutant ", colnames(pBAD_mat)[mutant])

	# construct AraC_data to fit
	ydat <- pBAD_mat[, mutant] # mutant, y-axis
	tdat <- pBAD_mat[, 1] # aTc, x-axis
	AraC_data <- data.frame(y=ydat, Lara=tdat)

	# define starting point and range of the parameters to feed into nls
	a1 = min(ydat) / 10
	a2 = min(ydat) * 2
	b1 = max(ydat) / 100
	b2 = max(ydat) * 2
	lower_bound <- c(a=a1, b=b1, k=0.001, n=1)
	upper_bound <- c(a=a2, b=b2, k=1, n=4)
	init_val <- c(a=a1, b=b1, k=0.1, n=2)
	obj_w <- 1 / ydat

	# use NLS (non-linear least square) to fit 4 parameters of each mutant with the experimental data
	anlsb1 <- try(nls(AraC_sigmoid, start=init_val, lower=lower_bound, upper=upper_bound, data=AraC_data, algorithm="port", weights=obj_w))
	AraC_data$predict <- predict(anlsb1, interval="predict") # append prediction

	message("AraC_data:")
	print(AraC_data)

	# extract model coefficients
	coefficients = coef(anlsb1)
	pBAD_parameter[, mutant - 1] = coefficients # populate the coefficients to pBAD_parameter

	# save predictions
	pBAD_pred[, mutant] = AraC_data$predict

	cat("\n")
}

# By now, all parameter values were found and stored in the matrix pBAD_parameter
# In the matrix pBAD_parameter, rows represent parameters and columns represent mutants.
message("pBAD_parameter:")
print(pBAD_parameter)

# save the predictions
write.table(pBAD_pred, file=file.path(project_data_dir, "pBAD_pred.txt"), row.names=FALSE, col.names=colnames(pBAD_mat), sep="\t")


#####################
# pBAD_pTET cascade #
#####################

k_convert_gfp_to_tetR = 100
pTET_mut_name_list = colnames(pTET_mat)
pBAD_mut_name_list = colnames(pBAD_mat)

# create empty zero matrix to populate the final parameters
final_mat <- matrix(0, nrow=nrow(pBAD_pTET_mat), ncol=6)
colnames(final_mat) <- c("10ng/ml aTc (exp)",
						 "0.1% ara + 10ng/ml aTc (exp)",
						 "100ng/mL aTc (exp)",
						 "10ng/ml aTc (sim)",
						 "0.1% ara + 10ng/ml aTc (sim)",
						 "100ng/mL aTc (sim)")

# for each cascade
for (cascade in (1:nrow(pBAD_pTET_mat))) {
	# pTET mutant
	pTET_mut = pBAD_pTET_mat[cascade, 2]
	pTET_idx = 0

	for (i in 2:length(pTET_mut_name_list)) {
		if (pTET_mut_name_list[i] == pTET_mut)
			pTET_idx = i - 1
	}

	if (pTET_idx == 0)
		stop(sprintf("%s not found", pTET_mut))

	# pBAD mutant
	pBAD_mut = pBAD_pTET_mat[cascade, 1]
	pBAD_idx = 0

	for (i in 2:length(pBAD_mut_name_list)) {
		if (pBAD_mut_name_list[i] == pBAD_mut)
			pBAD_idx = i - 1
	}

	if (pBAD_idx == 0)
		stop(sprintf("%s not found", pBAD_mut))

	# load previously found parameters
	a_pBAD = pBAD_parameter[1, pBAD_idx]
	b_pBAD = pBAD_parameter[2, pBAD_idx]
	k_pBAD = pBAD_parameter[3, pBAD_idx]
	n_pBAD = pBAD_parameter[4, pBAD_idx]

	a_pTET = pTET_parameter[1, pTET_idx]
	b_pTET = pTET_parameter[2, pTET_idx]
	k_pTET = pTET_parameter[3, pTET_idx]
	n_pTET = pTET_parameter[4, pTET_idx]

	######################
	# AraC_sigmoid model #
	######################
	# model to predict the TetR where there is no L-arabinose.
	# this model is the AraC_sigmoid model but we set Lara = 0,
	# and the output is normalized by a linear factor constant k_convert_gfp_to_tetR
	tetR_zero_ara = (a_pBAD + b_pBAD / (1 + (araC_const / k_pBAD) ^ n_pBAD)) / k_convert_gfp_to_tetR

	# model to predict the TetR where there 0.1% L-arabinose.
	# this model is the AraC_sigmoid model but we set Lara = 0.1,
	# and the output is normalized by a linear factor constant k_convert_gfp_to_tetR
	tetR_0_1_ara = (a_pBAD + b_pBAD / (1 + (araC_const / ((1 + (0.1 / k_Lara) ^ n_Lara) * k_pBAD)) ^ n_pBAD)) / k_convert_gfp_to_tetR

	######################
	# tetR_sigmoid model #
	######################
	# model to predict gfp when there is 10ng/ml aTc.
	# this model is the tetR_sigmoid model but we set aTc = 10,
	# and we set tetR_const = tetR_zero_ara as tetR is the output of the first module (i.e. AraC_sigmoid model)
	gfp_10ng_aTc = a_pTET + b_pTET / (1 + (tetR_zero_ara / ((1 + (10 / k_aTc) ^ n_aTc) * k_pTET)) ^ n_pTET)

	# model to predict the gfp when there is 10ng/ml aTc and 0.1% L-arabinose. (0.1% ara + 10ng/ml aTc)
	# this model is the tetR_sigmoid model but we set aTc = 10,
	# and we set tetR_const = tetR_0_1_ara as tetR is the output of the first module (i.e. AraC_sigmoid model) when there is 0.1% L-arabinose,
	gfp_0_1_Lara_10ng_aTc = a_pTET + b_pTET / (1 + (tetR_0_1_ara / ((1 + (10 / k_aTc) ^ n_aTc) * k_pTET)) ^ n_pTET)

	# model to predict gfp when there is 100ng/mL aTc.
	# this model is the tetR_sigmoid model but we set aTc = 100
	# and we set tetR_const = tetR_zero_ara
	gfp_100ng_aTc = a_pTET + b_pTET / (1 + (tetR_zero_ara / ((1 + (100 / k_aTc) ^ n_aTc) * k_pTET)) ^ n_pTET)

	# polulate final mat
	final_mat[cascade, ] = c(pBAD_pTET_mat[cascade, 3],
							 pBAD_pTET_mat[cascade, 4],
							 pBAD_pTET_mat[cascade, 5],
							 gfp_10ng_aTc,
							 gfp_0_1_Lara_10ng_aTc,
							 gfp_100ng_aTc)
}

message("final_mat:")
print(final_mat)

# save the final matrix
write.table(final_mat, file=file.path(project_data_dir, "final_mat.txt"), row.names=FALSE, sep="\t")

# find correlation
cor(final_mat[, 1], final_mat[, 4]) # 10ng/ml aTc
cor(final_mat[, 2], final_mat[, 5]) # 0.1% ara + 10ng/ml aTc
cor(final_mat[, 3], final_mat[, 6]) # 100ng/mL aTc
