\name{createParam}
\alias{createParamFiles}
\alias{writeParamFiles}
\alias{createParam}

\title{
createParam
}
\description{
The constructor function creates an \code{SYSargs2} S4 class object from
command-line string. Also, the function creates and saves the CWL param files.
The latest storages all the parameters required for running command-line software,
following the standard and specification defined on Common Workflow Language (CWL).
}
\usage{
createParamFiles(commandline, cwlVersion = "v1.1", class = "CommandLineTool",
                  results_path = "./results", module_load = "baseCommand",
                  file = "default", syntaxVersion = "v1",
                  writeParamFiles = TRUE, confirm = FALSE,
                  overwrite = FALSE, silent = FALSE)

writeParamFiles(sysargs, file = "default", overwrite = TRUE, silent = FALSE,
                syntaxVersion = "v1")
}

\arguments{
  \item{commandline}{
  string. Original command-line to create the CWL files from. Please see
  \code{Details} for more information.
}
  \item{cwlVersion}{
  string. The version of the Common Workflow Language. More information
  here: https://www.commonwl.org/.
}
  \item{class}{
  character. Name of Common Workflow Language Description class. The following
  is accepted: \code{CommandLineTool}.
}
  \item{results_path}{
  Path to the results folder. Default is \code{results}.
}
  \item{module_load}{
  string, Name of software to load by Environment Modules system. Default is
  "baseCommand", which creates a subfolder and two files: *.cwl and *.yml
  at \code{./param/cwl/}.
}
  \item{file}{
  character. Name and path of output files. If set to "default" then the name of
  the output files will have the pattern \code{<software>.cwl} and
  \code{<software>.yml}, where <software> will be what \code{baseCommand(x)}
  returns, when \code{x} is an object of class \code{SYSargs2}.
  Also, it creates a subfolder at \code{./param/cwl/} with name <software>.
}
  \item{syntaxVersion}{
  character. One of \code{"v1"} or \code{"v2"}, what CWL parsing syntax version to use.
  Default is to use the old version. The V2 version comes with more feature support,
  but has more syntax restrictions. See details.
}
  \item{writeParamFiles}{
  logical. If set to TRUE, it will write to file the content of the CWL
  files:\code{*.cwl} and \code{*.yml}. Default is \code{TRUE}.
}
  \item{confirm}{
  If set to \code{FALSE} and in an interactive section, it will prompt a question to
  proceed or not.
}
  \item{overwrite}{
   logical. If set to TRUE, existing files of the same name will be overwritten.
   Default is \code{FALSE}.
}
  \item{silent}{
  logical. If set to TRUE, all messages returned by the function will be
  suppressed. Default is \code{FALSE}.
}
  \item{sysargs}{
  Object of class \code{SYSargs2}. Output from the \code{createParamFiles}
  function.
}
}
\details{

    \subsection{Version 1 syntax}{

    - First line of the command-line object will be treated as the
    \code{baseCommand};

    - For argument lines (starting from the second line), any word before the first
      space with leading `-` or `--` in each will be treated as a prefix, like
      \code{-S} or \code{--min}. Any line without this first word will be treated
      as no prefix;

    - All defaults are placed inside \code{<...>};

    - First argument is the input argument type. \code{F} for "File", \code{int}
      for integer, \code{string} for string;

    - Optional: use the keyword \code{out} followed the type with a \code{,} comma
      separation to indicate if this argument is also a CWL output;

    - Then, use \code{:} to separate keywords and default values, any non-space
      value after the \code{:} will be treated as the default value;

    - If any argument has no default value, just a flag, like \code{--verbose},
      there no need to add any \code{<...>}.

    - The \code{\\} is not required, however for consistency it is recommended to add.

    }

    \subsection{Version 2 syntax}{

    - First line of the command-line object will be treated as the
    \code{baseCommand};

    - Each line specifies one argument and its default value.

    - Each line is composed with exact 2 \code{;} to seprate 3 columns.

    - Text before first \code{;} will be will used as prefix/names.
    If it starts with keyword “\code{p:}”, anything after “\code{p:}”
    and before the first \code{;} will be used as prefix, and the name
    of this position will be the prefix but with leading dash(s) ‘\code{-}’,
    ‘\code{–}’ removed. If there is any duplication, a number index will be
    added to the end. If there is no keyword “\code{p:}” before first \code{;},
    all text before first \code{;} will be the name.

    - If there is keyword “\code{p:}” before first \code{;} but nothing before
    and after the second \code{;}, this position will be treated as CWL argument
    instead of input.

    - Text between first and second \code{;} is type. Must be one of File, Directory,
    string, int, double, float, long, boolean.

    - Text after second \code{;} and before \code{\\} or end of the line is the
    default value. If it starts with keyword “\code{out}” or “\code{stdout}”, this
    position will also be added to outputs or standard output.

    - There is only 1 position with “\code{stdout}” allowed and usually it is the
    last position arguement.

    - Ending with “\code{\\}” is recommended but not required.

    }
}
\value{
\code{SYSargs2} object
}
\references{
For more details on \code{CWL}, please consult the following
page: https://www.commonwl.org/
}
\author{
Le Zhang and Daniela Cassol
}
\seealso{
\code{writeParamFiles}
\code{printParam}
\code{subsetParam}
\code{replaceParam}
\code{renameParam}
\code{appendParam}
\code{loadWorkflow}
\code{renderWF}
\code{showClass("SYSargs2")}
}

\examples{
## syntax version 1  example
command <- "
hisat2 \
    -S <F, out: ./results/M1A.sam> \
    -x <F: ./data/tair10.fasta> \
    -k <int: 1> \
    -min-intronlen <int: 30> \
    -max-intronlen <int: 3000> \
    -threads <int: 4> \
    -U <F: ./data/SRR446027_1.fastq.gz> \
    --verbose
"
cmd <- createParam(command, writeParamFiles=FALSE)
cmdlist(cmd)

## syntax version 2  example
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
cmdlist(cmd2)

}

\keyword{ utilities }
