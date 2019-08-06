#############################################################################
##  WeylModules Package
##  S. R. Doty
#############################################################################

SetPackageInfo( rec(
PackageName := "WeylModules",
Subtitle := "A GAP Package",
Version := "2.0",
Date := "5 August 2019",
ArchiveURL := 
          "http://doty.math.luc.edu/WeylModules",
ArchiveFormats := ".zip",

Persons := [
  rec( 
    LastName      := "Doty",
    FirstNames    := "S R",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "doty@math.luc.edu",
    WWWHome       := "http://doty.math.luc.edu",
    PostalAddress := Concatenation( [
                       "Department of Mathematics & Statistics\n",
                       "Loyola University Chicago\n",
                       "Chicago, IL 60626 USA" ] ),
    Place         := "Chicago",
    Institution   := "Loyola University Chicago"
  )
],

Status := "other",

README_URL := 
  "http://doty.math.luc.edu/WeylModules/README",
PackageInfoURL := 
  "http://doty.math.luc.edu/WeylModules/PackageInfo.g",


AbstractHTML := 
  "This package provides a collection of functions for computing with \
   Weyl modules for a semisimple simply connected algebraic group.",

PackageWWWHome := "http://doty.math.luc.edu/WeylModules",
               

PackageDoc := rec(
  # use same as in GAP            
  BookName  := "WeylModules",
  # format/extension can be one of .zoo, .tar.gz, .tar.bz2, -win.zip
  Archive := 
      "http://doty.math.luc.edu/WeylModules",
  ArchiveURLSubset := ["doc", "htm"],
  HTMLStart := "htm/chapters.htm",
  PDFFile   := "doc/manual.pdf",
  # the path to the .six file used by GAP's help system
  SixFile   := "doc/manual.six",
  # a longer title of the book, this together with the book name should
  # fit on a single text line (appears with the '?books' command in GAP)
  # LongTitle := "Elementary Divisors of Integer Matrices",
  LongTitle := "Weyl Modules Package",
  # Should this help book be autoloaded when GAP starts up? This should
  # usually be 'true', otherwise say 'false'. 
  Autoload  := true
),
Dependencies := rec( GAP:= ">=4.4", 
  NeededOtherPackages := [],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

AvailabilityTest := ReturnTrue,

BannerString := Concatenation( 
  "----------------------------------------------------------------\n",
  "Loading Weyl Module Package, Version ", ~.Version, "\n",
#  "by ", ~.Persons[1].FirstNames, " ", ~.Persons[1].LastName,
#        " (", ~.Persons[1].WWWHome, ")\n",
#  "   ", ~.Persons[2].FirstNames, " ", ~.Persons[2].LastName,
#        " (", ~.Persons[2].WWWHome, ")\n",
  "For help, type: ?WeylModule package \n",
  "----------------------------------------------------------------\n" ),

Autoload := false,

Keywords := ["Weyl modules", "algebraic groups"]
));


