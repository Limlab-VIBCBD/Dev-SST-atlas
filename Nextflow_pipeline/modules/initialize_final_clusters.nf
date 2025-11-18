process initialize {

    input:
      val outpath
      val tmppath


    script:
    """
    #!/usr/bin/env Rscript
    if(!dir.exists('${outpath}')){dir.create('${outpath}')}
    if(!dir.exists('${tmppath}')){dir.create('${tmppath}')}
    file.create("${tmppath}/remove_cells.txt")
    """
}