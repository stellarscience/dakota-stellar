\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{dakota-article}[03/14/2017 v0.1 Dakota modified article class]

% Passes and class options to the underlying article class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\newif\ifmy@dakotalogo
\DeclareOption{dakotalogo}{\my@dakotalogotrue}
\DeclareOption{nologo}{\my@dakotalogofalse}
\ExecuteOptions{nologo}
\ProcessOptions\relax

\ProcessOptions

% Load LaTeX's article class with the `twoside option
\LoadClass[twoside]{article}

% Redefine the page margins
% TODO: Adjust margins to your liking
\RequirePackage[left=3.0cm, right=3.0cm, top=3.0cm, bottom=3.0cm]{geometry}

% Remove the numbers from all the headings (\section, \subsection, etc.)
%\setcounter{secnumdepth}{-1}

% Define colors
\RequirePackage{xcolor}
\definecolor{darkblue}{RGB}{39, 92, 115}
\definecolor{lightblue}{RGB}{78,184,210}
\definecolor{midblue}{RGB}{54,125,155}

% Color section headings
\RequirePackage{sectsty}
\sectionfont{\color{midblue}}
\subsectionfont{\color{midblue}}

% Color caption and table labels, e.g. Figure 1:
\RequirePackage{caption}
\captionsetup[figure]{labelfont={color=midblue}}
\captionsetup[table]{labelfont={color=midblue}}

\newcommand*\samethanks[1][\value{footnote}]{\footnotemark[#1]}

\newcommand\@shortauthor{}
% define \theshortauthor to what is given
\newcommand\shortauthor[1]{\renewcommand\@shortauthor{#1}}


\RequirePackage{fancyhdr}% http://ctan.org/pkg/fancyhdr
\setlength{\headheight}{21pt}%
\fancypagestyle{myheader}{%
  \fancyhf{}% Clear all headers/footers
  \fancyhead[LO]{\raisebox{-.5\baselineskip}[0pt][0pt]{\bf \color{midblue} \@shorttitle}}
  \fancyhead[LE]{\raisebox{-.5\baselineskip}[0pt][0pt]{\bf \color{midblue} \thepage}}
  \fancyhead[RE]{\raisebox{-.5\baselineskip}[0pt][0pt]{\bf \color{midblue} \@shortauthor}}
  \fancyhead[RO]{\raisebox{-.5\baselineskip}[0pt][0pt]{\bf \color{midblue} \thepage}}
  \renewcommand{\headrulewidth}{2pt}% 2pt header rule
  \renewcommand{\headrule}{\hbox to\headwidth{%
    \color{gray}\leaders\hrule height \headrulewidth\hfill}}
  \renewcommand{\footrulewidth}{0pt}% No footer rule
}

% Define style of embedded code boxes
\RequirePackage[most]{tcolorbox}
\tcbuselibrary{listings,skins}

\lstdefinestyle{mycodestyle}{
language=C++,
%backgroundcolor=\color{gray!20!white},
basicstyle=\footnotesize,
keywordstyle=\color{blue},
                stringstyle=\color{red},
                commentstyle=\color{red},
                morecomment=[l][\color{magenta}]{\#},
}

\newtcblisting{codelisting}[2][]{
    arc=0pt, outer arc=0pt,
    listing only, 
    listing style=mycodestyle,
    title=\Large#2,
    colbacktitle=gray!60,
    title style={left color=darkblue, right color=lightblue!70, shading angle=45},
interior style={left color=gray!20,right color=white,shading angle=135},
enhanced, drop shadow,
    #1
    }

\newcommand\@shorttitle{}
\newcommand\shorttitle[1]{\renewcommand\@shorttitle{#1}}

\newcommand\@funding{}
\newcommand\funding[1]{\renewcommand\@funding{#1}}

\newcommand\@corauthor{}
\newcommand\corauthor[1]{\renewcommand\@corauthor{#1}}

\newcommand\@coremail{}
\newcommand\coremail[1]{\renewcommand\@coremail{#1}}

\newcommand\nomarkfootnote[1]{%
  \let\thefootnote\relax\footnotetext{#1}%
  \let\thefootnote\svthefootnote%
}

\def\@maketitle{%
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
  \let \footnote \thanks
{
\begin{tikzpicture}[remember picture,overlay]
\coordinate [below=2.5cm] (midpoint) at (current page.north);
\node [name=outercolourbar,
anchor=base,
rectangle, left color=darkblue, right color =darkblue!40,
text = white,
minimum width=\paperwidth,
minimum height=1.5cm,
text width=0.6\paperwidth] at (midpoint) {};
\node [name=colourbar,
anchor=base,
rectangle,left color=darkblue, right color=lightblue!70,
shading = axis,shading angle=45,
text = white,
minimum width=\paperwidth,
minimum height=1.3cm,
text width=0.6\paperwidth] at (midpoint) {};
\ifmy@dakotalogo
\node [right=7cm,text = white,] at  (colourbar.west) {\LARGE
  \@shorttitle \par};
\else
\node [text = white, minimum width=0.8\paperwidth,] at (midpoint) {\LARGE
  \@shorttitle \par};
\fi
% Define the point where the logo will go
\coordinate [right=4cm] (logo) at (colourbar.west);
% Set coordinate system origin
\begin{scope}[shift=(logo)]
  % Include the logo
\ifmy@dakotalogo
\node {\includegraphics[width=4cm]{DAKOTA_Arrow}};
\fi
\end{scope}
\end{tikzpicture}
}

    {\LARGE \color{midblue}\@title \par}%
    \vskip 1.5em%
    {\large
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        {\color{darkblue}\@author}
      \end{tabular}\par}%
    \vskip 1em%
    {\large \color{darkblue}\@date    \nomarkfootnote{{\color{darkblue}Disclosure}: \@funding}
\nomarkfootnote{{\color{darkblue}Correspondence author:} 
{\@corauthor \hskip 1em(\@coremail)}} }%
  \end{center}%
  \par
  \vskip 1.5em}

\usepackage{amsthm}
\renewenvironment{proof}{{\bfseries {\color{midblue}Proof}}}{\qed}
% declare a new theorem style
\newtheoremstyle{mystyle}%
{3pt}% Space above
{3pt}% Space below 
{\itshape}% Body font
{}% Indent amount
{\bfseries\color{midblue}}% Theorem head font
{.}% Punctuation after theorem head
{.5em}% Space after theorem head 
{}% Theorem head spec (can be left empty, meaning ‘normal’)
\theoremstyle{mystyle}

\usepackage{enumitem}
\setlist[itemize]{label=\textcolor{midblue}{\textbullet}}

\newcommand{\keywords}[1]{{{\small\bfseries{\color{darkblue}Keywords.}} #1}}
\newcommand{\amsid}[1]{{{\small\bfseries{\color{darkblue}AMS
      classifications.}} #1}}
\renewcommand*\abstractname{\color{darkblue} Abstract}

\endinput
