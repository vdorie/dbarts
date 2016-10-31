## based sorta on parallel::detectCores

guessNumCores <- function(logical = FALSE) {
  if (is.na(logical)) stop("logical cannot be NA")
  if (length(logical) != 1) stop("logical must be length 1")
  
  cGuess <- .Call(C_dbarts_guessNumCores)
  if (logical == FALSE && !is.na(cGuess[1])) return(cGuess[1])
  if (logical == TRUE  && !is.na(cGuess[2])) return(cGuess[2])

  systems <- list(freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
                  irix  = c("hinv | grep Processors | sed 's: .*::'", "hinv | grep '^Processor '| wc -l"),
                  solaris = if (logical) "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l" else
                                         "/bin/kstat -p -m cpu_info | grep :core_id | cut -f2 | uniq | wc -l")
  for (i in seq(systems)) {
    if (length(grep(paste0("^", names(systems)[i]), R.version$os)) > 0) {
      for (command in systems[i]) {
        commandResult <- gsub("^ +", "", system(command, TRUE)[1])
        if (length(grep("^[1-9]", commandResult))) return(as.integer(commandResult))
      }
    }
  }
  NA_integer_
}
