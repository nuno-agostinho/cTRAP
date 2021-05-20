# Minimize RAM usage when comparing against NCI60's matrix correlation (GE/DSen)
# Nuno Agostinho, 13 Jan 2021

# nci60Whole <- loadExpressionDrugSensitivityAssociation("NCI60")
# origin     <- nci60Whole

# Subset data to quickly experiment
nci60 <- origin
attrs <- attributes(nci60)
nci60 <- nci60[1:5, 1:5]
attributes(nci60) <- c(attributes(nci60),
                       attrs[!names(attrs) %in% names(attributes(nci60))])
attr(nci60, "rownames") <- rownames(nci60)
attr(nci60, "colnames") <- colnames(nci60)

# Write NCI60's GE vs DE association as HDF5 -----------------------------------
# writeExpressionDrugSensitivityCorHDF5(nci60)
writeExpressionDrugSensitivityCorHDF5(nci60Whole)

# Read NCI60's GE vs DE association as HDF5 ------------------------------------
obj <- readExpressionDrugSensitivityCorHDF5()
all(nci60Whole == obj, na.rm = TRUE)

# Check if loaded object is similar to original one ----------------------------
printNonIdenticalAttrs <- function(a, b) {
    ns <- names(attributes(a))
    for (i in ns) {
        isIdentical <- identical(attr(a, i), attr(b, i))
        if (!isIdentical) isIdentical <- all(attr(a, i) == attr(b, i), na.rm=T)
        if (!isIdentical) isIdentical <- i == "date"
        if (!isIdentical) print(i)
    }
}
printNonIdenticalAttrs(nci60Whole, obj)

obj <- readExpressionDrugSensitivityCorHDF5(rows=1:5, cols=1:5)
obj <- readExpressionDrugSensitivityCorHDF5(rows=1:2, cols=1:2)
obj <- readExpressionDrugSensitivityCorHDF5(rows=2:4, cols=3:5)
obj
