library(randomForest)

# Mode function for categorical/factor columns
mode_value <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

# Initialize missing values
initialize_missing <- function(X) {
  X_imp <- X
  for(col_name in colnames(X_imp)) {
    missing_idx <- which(is.na(X_imp[[col_name]]))
    if(length(missing_idx) == 0) next
    
    if(is.numeric(X_imp[[col_name]])) {
      X_imp[missing_idx, col_name] <- quantile(
        X_imp[[col_name]], probs = 0.5, type = 0, na.rm = TRUE
      )
    } else {
      X_imp[missing_idx, col_name] <- mode_value(X_imp[[col_name]])
    }
  }
  return(X_imp)
}

# Iterative training imputation, return both imputed data + final RFs
impute_train_and_save_rfs <- function(X, maxiter = 10, tol = 1e-3) {
  X_imp <- initialize_missing(X)
  final_rf_list <- list()
  
  for(iter in 1:maxiter) {
    X_old <- X_imp
    rf_list_iter <- list()  # temporary RFs for this iteration
    
    for(col_name in colnames(X_imp)) {
      missing_idx <- which(is.na(X[[col_name]]))
      if(length(missing_idx) == 0) next
      
      predictor_cols <- setdiff(colnames(X_imp), col_name)
      
      # Refit RF on current imputed dataset
      rf_model <- randomForest(
        x = X_imp[, predictor_cols],
        y = X_imp[[col_name]],
        na.action = na.omit
      )
      rf_list_iter[[col_name]] <- rf_model
      
      # Predict missing entries
      X_imp[missing_idx, col_name] <- predict(rf_model, X_imp[missing_idx, predictor_cols])
    }
    
    # Update final RFs if convergence reached
    final_rf_list <- rf_list_iter
    
    # Convergence check
    diff <- sum((as.matrix(X_imp) - as.matrix(X_old))^2, na.rm = TRUE) /
      sum((as.matrix(X_old))^2, na.rm = TRUE)
    if(diff < tol) break
  }
  
  return(list(imputed_data = X_imp, final_rfs = final_rf_list))
}


impute_test_with_rfs <- function(rf_list, 
                                 X_test, 
                                 maxiter = 10, 
                                 tol = 1e-3) {
  X_imp <- initialize_missing(X_test)
  
  for(iter in 1:maxiter) {
    X_old <- X_imp
    
    for(col_name in names(rf_list)) {
      missing_idx <- which(is.na(X_test[[col_name]]))
      if(length(missing_idx) == 0) next
      
      predictor_cols <- setdiff(colnames(X_imp), col_name)
      rf_model <- rf_list[[col_name]]  # frozen RF from training
      
      # Predict missing entries
      X_imp[missing_idx, col_name] <- predict(rf_model, X_imp[missing_idx, predictor_cols])
    }
    
    # Convergence check
    diff <- sum((as.matrix(X_imp) - as.matrix(X_old))^2, na.rm = TRUE) /
      sum((as.matrix(X_old))^2, na.rm = TRUE)
    if(diff < tol) break
  }
  
  return(X_imp)
}

# Example
# impute_train <- impute_train_and_save_rfs(train)
# x_tr <-impute_train$imputed_data
# x_te <- impute_test_with_rfs(train_result$final_rfs, test)

