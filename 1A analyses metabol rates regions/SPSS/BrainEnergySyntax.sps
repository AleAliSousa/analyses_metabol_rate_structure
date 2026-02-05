* Encoding: UTF-8.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc slope MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: slope=col(source(s), name("slope"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("slope"))
  ELEMENT: point(position(rCMRGlc*slope))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc lambd MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: lambd=col(source(s), name("lambd"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("lambd"))
  ELEMENT: point(position(rCMRGlc*lambd))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc OPdiff MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: OPdiff=col(source(s), name("OPdiff"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("OPdiff"))
  ELEMENT: point(position(rCMRGlc*OPdiff))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc ageBeta MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: ageBeta=col(source(s), name("ageBeta"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("ageBeta"))
  ELEMENT: point(position(rCMRGlc*ageBeta))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc ageChange MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: ageChange=col(source(s), name("ageChange"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("ageChange"))
  ELEMENT: point(position(rCMRGlc*ageChange))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc devChange MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: devChange=col(source(s), name("devChange"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("devChange"))
  ELEMENT: point(position(rCMRGlc*devChange))
END GPL.

* Log


* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=logrCMRGlc slope MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: logrCMRGlc=col(source(s), name("logrCMRGlc"))
  DATA: slope=col(source(s), name("slope"))
  GUIDE: axis(dim(1), label("logrCMRGlc"))
  GUIDE: axis(dim(2), label("slope"))
  ELEMENT: point(position(logrCMRGlc*slope))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=logrCMRGlc lambd MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: logrCMRGlc=col(source(s), name("logrCMRGlc"))
  DATA: lambd=col(source(s), name("lambd"))
  GUIDE: axis(dim(1), label("logrCMRGlc"))
  GUIDE: axis(dim(2), label("lambd"))
  ELEMENT: point(position(logrCMRGlc*lambd))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=logrCMRGlc OPdiff MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: logrCMRGlc=col(source(s), name("logrCMRGlc"))
  DATA: OPdiff=col(source(s), name("OPdiff"))
  GUIDE: axis(dim(1), label("logrCMRGlc"))
  GUIDE: axis(dim(2), label("OPdiff"))
  ELEMENT: point(position(logrCMRGlc*OPdiff))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=logrCMRGlc ageBeta MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: logrCMRGlc=col(source(s), name("logrCMRGlc"))
  DATA: ageBeta=col(source(s), name("ageBeta"))
  GUIDE: axis(dim(1), label("logrCMRGlc"))
  GUIDE: axis(dim(2), label("ageBeta"))
  ELEMENT: point(position(logrCMRGlc*ageBeta))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=logrCMRGlc ageChange MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: logrCMRGlc=col(source(s), name("logrCMRGlc"))
  DATA: ageChange=col(source(s), name("ageChange"))
  GUIDE: axis(dim(1), label("logrCMRGlc"))
  GUIDE: axis(dim(2), label("ageChange"))
  ELEMENT: point(position(logrCMRGlc*ageChange))
END GPL.

* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=logrCMRGlc devChange MISSING=LISTWISE REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: logrCMRGlc=col(source(s), name("logrCMRGlc"))
  DATA: devChange=col(source(s), name("devChange"))
  GUIDE: axis(dim(1), label("logrCMRGlc"))
  GUIDE: axis(dim(2), label("devChange"))
  ELEMENT: point(position(logrCMRGlc*devChange))
END GPL.



DATASET ACTIVATE DataSet1.
* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc AstroDensity MISSING=LISTWISE REPORTMISSING=NO    
  /GRAPHSPEC SOURCE=INLINE
  /FITLINE TOTAL=YES.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: AstroDensity=col(source(s), name("AstroDensity"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("AstroDensity"))
  GUIDE: text.title(label("Simple Scatter with Fit Line of AstroDensity by rCMRGlc"))
  ELEMENT: point(position(rCMRGlc*AstroDensity))
END GPL.



* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc MicroDensity MISSING=LISTWISE REPORTMISSING=NO    
  /GRAPHSPEC SOURCE=INLINE
  /FITLINE TOTAL=YES.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: MicroDensity=col(source(s), name("MicroDensity"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("MicroDensity"))
  GUIDE: text.title(label("Simple Scatter with Fit Line of MicroDensity by rCMRGlc"))
  ELEMENT: point(position(rCMRGlc*MicroDensity))
END GPL.




* Chart Builder.
GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=rCMRGlc GliaNeurDensityRatio MISSING=LISTWISE 
    REPORTMISSING=NO
  /GRAPHSPEC SOURCE=INLINE
  /FITLINE TOTAL=YES.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: rCMRGlc=col(source(s), name("rCMRGlc"))
  DATA: GliaNeurDensityRatio=col(source(s), name("GliaNeurDensityRatio"))
  GUIDE: axis(dim(1), label("rCMRGlc"))
  GUIDE: axis(dim(2), label("GliaNeurDensityRatio"))
  GUIDE: text.title(label("Simple Scatter with Fit Line of GliaNeurDensityRatio by rCMRGlc"))
  ELEMENT: point(position(rCMRGlc*GliaNeurDensityRatio))
END GPL.
