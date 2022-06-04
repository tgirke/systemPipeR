\name{printParam2}
\alias{printParam2}
\alias{removeParam2}
\alias{replaceParam2}
\alias{renameParam2}
\alias{appendParam2}

\title{
Accessories function to modify the Command-line Version 2
}
\description{
Accessories function to modify the Command-line Version 2
}
\usage{
printParam2(sysargs, base = FALSE, args = FALSE, inputs = FALSE,
        outputs = FALSE, stdout = FALSE, raw_cmd = FALSE, all = TRUE)
appendParam2(sysargs, x, position = c("inputs", "args", "outputs"),
        after = NULL, verbose = FALSE)
replaceParam2(sysargs, x, index=NULL,
              position = c("inputs", "baseCommand", "args", "outputs", "stdout"),
              verbose = FALSE)
removeParam2(sysargs, index=NULL, position = c("inputs", "args", "outputs", "stdout"),
             verbose = FALSE)
renameParam2(sysargs, index=NULL, new_names,
             position = c("inputs", "args", "outputs", "stdout"), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sysargs}{
  Object of class \code{SYSargs2}. Output from the \code{createParamFiles} function.
}
  \item{base}{
  logical, print out base command information?
}
  \item{args}{
  logical, print out arguments information?
}
  \item{inputs}{
  logical, print out inputs information?
}
  \item{outputs}{
  logical, print out outputs information?
}
  \item{stdout}{
  logical, print out stdout information?
}
  \item{raw_cmd}{
  logical, print out parsed raw command information?
}
  \item{all}{
  logical, print out all base command, arguments, inputs,
  outputs, and raw command information?
  Turn this to \code{TRUE} will overwrite all above to \code{TRUE}.
  If you need to print out selected positions, turn this to \code{FALSE} and turn
  other positions to \code{TRUE}.
}
  \item{position}{
  string, one of the positions to apply a modification.
  For \code{appendParam2}: "inputs", "args", "outputs",
  for \code{replaceParam2}: "inputs", "baseCommand", "args", "outputs", "stdout",
  for \code{removeParam2}, \code{renameParam2}: "inputs", "args", "outputs", "stdout".
}
  \item{index}{
  numeric or character vector, index of items to remove or rename item(s)
  inside the position you choose, in \code{removeParam2},  \code{renameParam2},
  \code{replaceParam2}.
}
  \item{after}{
  integer, in \code{appendParam2}, after which current item you want to append
  the new items in the position you have selected? For example, if you want
  the new item to be the first, then use \code{0}, if last, then use
  \code{999} to make sure it goes to the last one.
}
  \item{x}{
  named list or string, new items to replace or append in different positions. Replace list length
  must be the same as \code{index}. Different positions will have different requirements.
  See details of the \code{x} format requirements.
}
  \item{new_names}{
  character vector, new names that you wish to replace the old names. \code{new_names} vector length
  must be the same as \code{index} in \code{renameParam2}.
}
  \item{verbose}{
  logical, show addtional information during/after operation? for example, print the new changes.
}
}

\details{
- \code{printParam2}: prints its arguments defined by \code{position} and \code{index}.

- \code{removeParam2}: removes items in certain positions you select.

- \code{replaceParam2}: replaces the values in command-line with indices given in list by those given in values

- \code{renameParam2}: rename the names of items in certain position.

- \code{appendParam2}: Add arguments to the original command line. Adding new basecommand or
  standard out is not allowed.

  \subsection{x format}{
    - If \code{x} is a character, it requires exact 3 semi-colons \code{;} to separate the string
    in to 4 columns. Values before the third column is the same as \code{createParam} inputs,
    first column: prefix/argument name, second column: type, third column: default value.
    The fourth column (new): numeric, index of the new item, this will be translated into
    \code{position} entries in CWL.

    - If \code{x} is a list, it must be named. Following items must be included in list:
     \code{preF}, \code{type}, \code{value}, \code{index}. They refer to prefix, param type,
     default value, and position index correspondingly.
  }
}
\value{
\code{SYSargs2} object
}
\references{
For more details on \code{CWL}, please consult the following page: https://www.commonwl.org/
}
\author{
Le Zhang and Daniela Cassol
}
\seealso{
\code{writeParamFiles}
\code{createParamFiles}
\code{loadWorkflow}
\code{renderWF}
\code{showClass("SYSargs2")}
}

\examples{
command2 <- '
mycmd2 \
    p: -s; File; sample1.txt \
    p: -s; File; sample2.txt \
    p: --c; ; \
    p: -o; File; out: myout.txt \
    ref_genome; File; a.fasta \
    p: --nn; int; 12 \
    mystdout; File; stdout: abc.txt
'
cmd2 <- createParam(command2, syntaxVersion = "v2", writeParamFiles=FALSE)
# string format
new_cmd <- 'p: -abc; string; abc; 7'
cmd2 <- appendParam2(cmd2, x = new_cmd, position = "inputs")
printParam2(cmd2, all = FALSE, inputs = TRUE, raw_cmd = TRUE)
# list format
new_cmd <- list(name = "new_arg", preF = "--foo", index = "8")
cmd2 <- appendParam2(cmd2, x = new_cmd, position = "args")
printParam2(cmd2, all = FALSE, args = TRUE, raw_cmd = TRUE)
# rename
cmd2 <- renameParam2(cmd2, "new_name_arg", index = "new_arg", position = "args")
printParam2(cmd2, all = FALSE, args = TRUE, raw_cmd = TRUE)
# remove
cmd2 <- removeParam2(cmd2, index = "new_name_arg", position = "args")
printParam2(cmd2, all = FALSE, args = TRUE, raw_cmd = TRUE)
}

\keyword{ utilities }

