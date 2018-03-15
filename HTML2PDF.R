basename <- "JDD"
decktape <- "../../useR-2017/decktape-1.0.0/"

system(glue::glue(
  "{decktape}phantomjs {decktape}decktape.js", 
  " {basename}.html {basename}.pdf", 
  " --load-pause 1000", 
  " -s 1504x1129"
))
