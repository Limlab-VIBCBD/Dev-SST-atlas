process copy_outputs {

    input:
      val outpath
      val tmppath
      path clusters
      path obj_with_clustering
      path subsample_ident_file
      path fullsample_ident_file
      path nPCs

   script:
    """
    cp ${clusters} ${tmppath}/
    cp ${obj_with_clustering} ${tmppath}/
    cp ${subsample_ident_file} ${tmppath}/
    cp ${fullsample_ident_file} ${tmppath}/
    cp ${nPCs} ${tmppath}/
    """
}