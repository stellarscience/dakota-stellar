file(READ "${infile}" contents)
string(REGEX REPLACE ${find} ${replace} contents "${contents}")
file(WRITE "${outfile}" "${contents}")
file(REMOVE "${infile}")
