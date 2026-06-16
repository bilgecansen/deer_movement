rc <- tibble::tribble(
  ~from , ~to  , ~id , ~label                ,
   4000 , 4301 ,   1 , "forest"              ,
   6400 , 6451 ,   2 , "forested wetlands"   ,
   6300 , 6331 ,   3 , "lowland scrub/scrub"
)
#Permissive: rc only needs to cover what you care about; everything else is NA
r <- assign_landcover(src = wiscland2, rc = rc, na_unmapped = TRUE)
