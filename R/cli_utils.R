# -i INPUT, --input INPUT
get_arg <- function(args, flags, default = NULL, required = FALSE) {
    i <- which(args %in% flags)

    if (length(i) == 0 || i[1] == length(args) || grepl("^--?", args[i[1] + 1])) {
        if (required) {
            stop(sprintf("Missing required argument: %s", paste(flags, collapse = " or ")))
        } else {
            return(default)
        }
    }

    args[i[1] + 1]
}

# --input 0.5
get_numeric_arg <- function(args, flags, default = NULL) {
    val <- get_arg(args, flags)
    if (is.null(val)) {
        return(default)
    }
    as.numeric(val)
}

# --input 7
get_integer_arg <- function(args, flags, default = NULL) {
    val <- get_arg(args, flags)
    if (is.null(val)) {
        return(default)
    }
    as.integer(val)
}

# --input a,b,c
get_list_arg <- function(args, flag) {
    val <- get_arg(args, flag)
    if (is.null(val)) {
        return(NULL)
    }
    strsplit(val, ",")[[1]]
}

# --input 1,2,3
get_numeric_list_arg <- function(args, flag) {
    val_list <- get_list_arg(args, flag)
    if (is.null(val_list)) {
        return(NULL)
    }
    as.numeric(val_list)
}

# --verbose
store_true <- function(args, flags) {
    any(flags %in% args)
}

# --quiet
store_false <- function(args, flags) {
    !any(flags %in% args)
}

# --fixed_column abc --> "abc" (str); --fixed_column 123 --> 123 (int); --fixed_column "123" --> 123 (int); --fixed_column str:123 --> "123" (str)
get_fixed_column <- function(args, flag, default = NULL) {
    val <- get_arg(args, flag)
    if (is.null(val)) {
        return(default)
    }

    if (startsWith(val, "str:")) {
        return(sub("^str:", "", val))
    } else if (grepl("^\\d+$", val)) {
        return(as.integer(val))
    } else {
        return(val)
    }
}

# --inputs A B C --outputs ... (vector) OR --inputs A=1 B=2 C=3 outputs ... (named vector)
get_multi_arg <- function(args, flags, required = FALSE) {
    i <- which(args %in% flags)
    
    if (length(i) == 0) {
        if (required) stop(sprintf("Missing required argument: one of %s", paste(flags, collapse = ", ")))
        return(NULL)
    }
    
    i <- i[1]
    
    # Collect values until next flag
    vals <- character()
    j <- i + 1
    while (j <= length(args) && !grepl("^-", args[j])) {
        vals <- c(vals, args[j])
        j <- j + 1
    }
    
    if (required && length(vals) == 0) {
        stop(sprintf("No values provided for required argument: %s", flags))
    }
    
    # Parse into named/unnamed values
    keys <- character(length(vals))
    vals_out <- character(length(vals))
    
    for (k in seq_along(vals)) {
        if (grepl("=", vals[k])) {
            kv <- strsplit(vals[k], "=", fixed = TRUE)[[1]]
            keys[k] <- kv[1]
            vals_out[k] <- kv[2]
        } else {
            vals_out[k] <- vals[k]
        }
    }
    
    # Apply names only if at least one key is non-empty
    if (any(nzchar(keys))) {
        names(vals_out) <- keys
    }
    
    return(vals_out)
}

