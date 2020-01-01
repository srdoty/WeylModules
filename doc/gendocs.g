path:=Directory(".");; # current directory contains the xmal file(s)
main:="weyl.xml";; # main xml file
files := ["../lib/weylmod.gd", "../lib/weylmod.gi"];; # should not be needed
bookname := "Weyl";;

MakeGAPDocDoc(path, main, files, bookname,"MathJax");;

# Run the following command-line command:
# $ gap gendocs.g
# from this folder (/doc) to generate all the documentation.
