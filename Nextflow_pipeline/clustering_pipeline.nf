#!/usr/bin/env nextflow

//nextflow.preview.recursion = true
//nextflow.enable.dsl=2

// Load modules
include {initialize} from './modules/initialize_final_clusters.nf'
include {integration} from './modules/integration.nf'
include {prepare_PC_5} from './modules/prepare_PC5.nf'
include {compute_PC_5} from './modules/compute_PC5.nf'
include {prepare_PC_1} from './modules/prepare_PC1.nf'
include {compute_nPCs_and_clusterings} from './modules/compute_nPCs_and_clusterings.nf'
include {resolutions_to_test} from './modules/resolutions_to_test.nf'
include {subsampling_and_clustering} from './modules/subsampling_and_clustering.nf'
include {identify_stable_clusters} from './modules/identify_stable_clusters.nf'
include {compute_final_clusters} from './modules/compute_final_clusters.nf'
include {copy_outputs} from './modules/copy_outputs.nf'
include {update_cluster_id} from './modules/update_cluster_id.nf'
include {split_object} from './modules/split_object.nf'
include {integration_iteration} from './modules/integration_iteration.nf'
include {prepare_PC_5_iteration} from './modules/prepare_PC5_iteration.nf'
include {compute_PC_5_iteration} from './modules/compute_PC5_iteration.nf'
include {prepare_PC_1_iteration} from './modules/prepare_PC1_iteration.nf'
include {compute_nPCs_and_clusterings_iteration} from './modules/compute_nPCs_and_clusterings_iteration.nf'
include {resolutions_to_test_iteration} from './modules/resolutions_to_test_iteration.nf'
include {subsampling_and_clustering_iteration} from './modules/subsampling_and_clustering_iteration.nf'
include {identify_stable_clusters_iteration} from './modules/identify_stable_clusters_iteration.nf'
include {compute_final_clusters_iteration} from './modules/compute_final_clusters_iteration.nf'
include {reassign_discarded_cells} from './modules/reassign_discarded_cells.nf'

// Clustering pipeline workflow
workflow run_pipeline {
   take:
    object

   main:
    integration(object,params.tmp_path,params.reproducibility)
      k = channel.of(1..10)
    prepare_PC_5(integration.out.data_use, integration.out.dgem_kfold, k)
    compute_PC_5(integration.out.data_use,prepare_PC_5.out.collect(),params.reproducibility)
    prepare_PC_1(integration.out.data_use,compute_PC_5.out.nPCs_5, compute_PC_5.out.dgem_kfold, k)
    compute_nPCs_and_clusterings(integration.out.obj_with_pca,compute_PC_5.out.nPCs_5,prepare_PC_1.out.collect(),params.min_resolution,params.max_resolution)
    resolutions_to_test(compute_nPCs_and_clusterings.out.obj_with_clustering,params.min_resolution,params.max_resolution,params.tmp_path)
        res = resolutions_to_test.out.resolutions.splitCsv().map { line -> line[0] }
        sub = channel.of(1..params.n_subsampling)
        combined_ch = res.combine(sub)
    subsampling_and_clustering(compute_nPCs_and_clusterings.out.obj_with_clustering,compute_nPCs_and_clusterings.out.nPCs,combined_ch,params.perc_sub,params.reproducibility)
    identify_stable_clusters(subsampling_and_clustering.out.collect(), compute_nPCs_and_clusterings.out.obj_with_clustering,params.jaccard_cutoff,params.percent_cutoff)
    compute_final_clusters(compute_nPCs_and_clusterings.out.obj_with_clustering,compute_nPCs_and_clusterings.out.nPCs,identify_stable_clusters.out.stable_resolution,params.minSize,params.DEscore_cutoff)
    copy_outputs(params.out_path,params.tmp_path,compute_final_clusters.out.clusters, compute_nPCs_and_clusterings.out.obj_with_clustering,identify_stable_clusters.out.subsample_ident_file,identify_stable_clusters.out.fullsample_ident_file,compute_nPCs_and_clusterings.out.nPCs)
    compute_final_clusters.out.num_clusters.map{ file -> file.text.trim() }.set{ num_cl }
    compute_nPCs_and_clusterings.out.obj_with_clustering.view { it -> "Pipeline produced file $it" }
    num_cl.view {it -> "Pipeline produced $it clusters" }
    
   emit:
   clusters = compute_final_clusters.out.clusters
   object = compute_nPCs_and_clusterings.out.obj_with_clustering
   num_clusters = num_cl
}

// Clustering pipeline for iterations workflow
workflow run_pipeline_iteration {
   take:
    object_cl

   main:
    integration_iteration(object_cl,params.tmp_path,params.reproducibility)
    data_use_ch = integration_iteration.out.map{ t -> tuple(t[1],t[2],t[3])}.collect(flat: false).flatMap() // this is done to wait integration_iteration to finish runnning on the whole channel befor running next process
    k = channel.of(1..10)
    data_k_combined = data_use_ch.combine(k)
    prepare_PC_5_iteration(data_k_combined)
    P5_ch = prepare_PC_5_iteration.out.collect(flat: false).flatMap().groupTuple(by: [1])
    data_use_ch = integration_iteration.out.map{ t -> tuple(t[1],t[3])}
    P5_ch = P5_ch.combine(data_use_ch, by: 1)
    compute_PC_5_iteration(P5_ch,params.reproducibility)
    data_k_combined_1_tmp = compute_PC_5_iteration.out.collect(flat: false).flatMap()
    data_k_combined_1 = data_k_combined_1_tmp.combine(k)
    prepare_PC_1_iteration(data_k_combined_1)
    P1_ch = prepare_PC_1_iteration.out.collect(flat: false).flatMap().groupTuple(by: [1])
    data_obj_ch = integration_iteration.out.map{ t -> tuple(t[0],t[3])}
    P1_ch = P1_ch.combine(data_obj_ch, by: 1)
    data_obj_ch = compute_PC_5_iteration.out.map{ t -> tuple(t[3],t[0],t[1])}
    P1_ch = P1_ch.combine(data_obj_ch, by: 0)
    compute_nPCs_and_clusterings_iteration(P1_ch,params.min_resolution,params.max_resolution)
    resolutions_to_test_iteration(compute_nPCs_and_clusterings_iteration.out, params.min_resolution,params.max_resolution)
    ch_expanded = resolutions_to_test_iteration.out.flatMap{object, nPC, cl, f -> def res = f.text.readLines()
                        res.collect{line -> tuple(line, object, nPC, cl)}}
    sub = channel.of(1..params.n_subsampling)
    combined_ch = ch_expanded.combine(sub)
    subsampling_and_clustering_iteration(combined_ch,params.perc_sub,params.reproducibility)
    subsampling_ch = subsampling_and_clustering_iteration.out.collect(flat: false).flatMap().groupTuple(by: 0)
    subsampling_ch=subsampling_ch.combine(compute_nPCs_and_clusterings_iteration.out, by: 0)
    identify_stable_clusters_iteration(subsampling_ch,params.jaccard_cutoff,params.percent_cutoff)
    compute_final_clusters_iteration(identify_stable_clusters_iteration.out.map{ t -> tuple(t[0],t[1],t[2],t[3])},params.minSize,params.DEscore_cutoff)
    final_objects = compute_final_clusters_iteration.out.flatMap{ obj, clusters, cl, file -> def num_cl=file.text.trim() 
                                                      num_cl.collect{line -> tuple(obj, clusters, cl, line)}}
    final_objects.view()
    final_objects.map{t -> t[3]}.view { it -> "Pipeline produced $it clusters" }

   emit:
    clusters = compute_final_clusters_iteration.out.map{t -> t[1]}.collect()
    object = compute_final_clusters_iteration.out.map{t -> t[0]}.collect()
    num_clusters = compute_final_clusters_iteration.out.map{t -> t[3]}.map{ file -> file.text.trim() }.collect()
}

// Run second iteration workflow
workflow run_iteration{
    take:
     out_ch

     main:
      input_split=out_ch.map{t -> tuple(t[0],t[1],t[2].toInteger())}.flatMap{obj,cluster,num -> def nn=(1..num).toList()
                                             nn.collect{n -> tuple(obj,cluster,n)}}
      update_cluster_id(input_split)
      input_split_updated = update_cluster_id.out.map{ object, clusters, cl, f ->
              def cl_updated = f.text.trim()
              tuple(object, clusters, cl, cl_updated)
          }
      split_object(input_split_updated,params.minSize)
      splitted_obj = split_object.out.collect(flat: false).flatMap()
      splitted_obj.view()
      splitted_obj = splitted_obj.map{ object, clusters,  f ->
              def size = f.text.trim()
              tuple(object, clusters, size)
          }
      splitted_obj=splitted_obj.filter{it -> it[2].toInteger() >= 2*(params.minSize+1)} 
      splitted_obj=splitted_obj.map{t -> tuple(t[0],t[1])}
      splitted_obj.view()
      run_pipeline_iteration(splitted_obj)
      
      c1=run_pipeline_iteration.out.object.flatMap { list -> list } 
      c2=run_pipeline_iteration.out.clusters.flatMap { list -> list } 
      c3=run_pipeline_iteration.out.num_clusters.flatMap { list -> list } 
      ch=c1.merge(c2).merge(c3) 
      out_ch=ch.filter{it -> it[2].toInteger() > 1} 

    emit:
      out_ch=out_ch
}

// Run third iteration workflow
workflow run_iteration_recurse{
    take:
     out_ch

     main:
      input_split=out_ch.map{t -> tuple(t[0],t[1],t[2].toInteger())}.flatMap{obj,cluster,num -> def nn=(1..num).toList()
                                             nn.collect{n -> tuple(obj,cluster,n)}}
      update_cluster_id(input_split)
      update_out=update_cluster_id.out.collect(flat: false).flatMap()
      input_split_updated = update_out.map{ object, clusters, cl, f ->
              def cl_updated = f.text.trim()
              tuple(object, clusters, cl, cl_updated)
          }
      split_object(input_split_updated,params.minSize) 
      splitted_obj = split_object.out.collect(flat: false).flatMap()
      splitted_obj.view()
      splitted_obj = splitted_obj.map{ object, clusters,  f ->
              def size = f.text.trim()
              tuple(object, clusters, size)
          }
      splitted_obj=splitted_obj.filter{it -> it[2].toInteger() >= 2*(params.minSize+1)} 
      splitted_obj=splitted_obj.map{t -> tuple(t[0],t[1])}
      splitted_obj.view()

      run_pipeline_iteration(splitted_obj)

      c1=run_pipeline_iteration.out.object.flatMap { list -> list } 
      c2=run_pipeline_iteration.out.clusters.flatMap { list -> list } 
      c3=run_pipeline_iteration.out.num_clusters.flatMap { list -> list } 
      ch=c1.merge(c2).merge(c3) 
      out_ch=ch.filter{it -> it[2].toInteger() > 1} 

    emit:
      out_ch=out_ch
}

// Run fourth iteration workflow
workflow run_iteration_recurse1{ 
    take:
     out_ch

     main:
      input_split=out_ch.map{t -> tuple(t[0],t[1],t[2].toInteger())}.flatMap{obj,cluster,num -> def nn=(1..num).toList()
                                             nn.collect{n -> tuple(obj,cluster,n)}}
      update_cluster_id(input_split)
      update_out=update_cluster_id.out.collect(flat: false).flatMap()
      input_split_updated = update_out.map{ object, clusters, cl, f ->
              def cl_updated = f.text.trim()
              tuple(object, clusters, cl, cl_updated)
          }
      split_object(input_split_updated,params.minSize) 
      splitted_obj = split_object.out.collect(flat: false).flatMap()
      splitted_obj.view()
      splitted_obj = splitted_obj.map{ object, clusters,  f ->
              def size = f.text.trim()
              tuple(object, clusters, size)
          }
      splitted_obj=splitted_obj.filter{it -> it[2].toInteger() >= 2*(params.minSize+1)} 
      splitted_obj=splitted_obj.map{t -> tuple(t[0],t[1])}
      splitted_obj.view()

      run_pipeline_iteration(splitted_obj)

      c1=run_pipeline_iteration.out.object.flatMap { list -> list } 
      c2=run_pipeline_iteration.out.clusters.flatMap { list -> list } 
      c3=run_pipeline_iteration.out.num_clusters.flatMap { list -> list } 
      ch=c1.merge(c2).merge(c3) 
      out_ch=ch.filter{it -> it[2].toInteger() > 1} 

    emit:
      out_ch=out_ch
}

// Run fifth iteration workflow
workflow run_iteration_recurse2{ 
    take:
     out_ch

     main:
      input_split=out_ch.map{t -> tuple(t[0],t[1],t[2].toInteger())}.flatMap{obj,cluster,num -> def nn=(1..num).toList()
                                             nn.collect{n -> tuple(obj,cluster,n)}}
      update_cluster_id(input_split)
      update_out=update_cluster_id.out.collect(flat: false).flatMap()
      input_split_updated = update_out.map{ object, clusters, cl, f ->
              def cl_updated = f.text.trim()
              tuple(object, clusters, cl, cl_updated)
          }
      split_object(input_split_updated,params.minSize) 
      splitted_obj = split_object.out.collect(flat: false).flatMap()
      splitted_obj.view()
      splitted_obj = splitted_obj.map{ object, clusters,  f ->
              def size = f.text.trim()
              tuple(object, clusters, size)
          }
      splitted_obj=splitted_obj.filter{it -> it[2].toInteger() >= 2*(params.minSize+1)} 
      splitted_obj=splitted_obj.map{t -> tuple(t[0],t[1])}
      splitted_obj.view()

      run_pipeline_iteration(splitted_obj)

      c1=run_pipeline_iteration.out.object.flatMap { list -> list } 
      c2=run_pipeline_iteration.out.clusters.flatMap { list -> list } 
      c3=run_pipeline_iteration.out.num_clusters.flatMap { list -> list } 
      ch=c1.merge(c2).merge(c3) 
      out_ch=ch.filter{it -> it[2].toInteger() > 1} 

    emit:
      out_ch=out_ch
}

// Run sixth iteration workflow
workflow run_iteration_recurse3{ 
    take:
     out_ch

     main:
      input_split=out_ch.map{t -> tuple(t[0],t[1],t[2].toInteger())}.flatMap{obj,cluster,num -> def nn=(1..num).toList()
                                             nn.collect{n -> tuple(obj,cluster,n)}}
      update_cluster_id(input_split)
      update_out=update_cluster_id.out.collect(flat: false).flatMap()
      input_split_updated = update_out.map{ object, clusters, cl, f ->
              def cl_updated = f.text.trim()
              tuple(object, clusters, cl, cl_updated)
          }
      split_object(input_split_updated,params.minSize) 
      splitted_obj = split_object.out.collect(flat: false).flatMap()
      splitted_obj.view()
      splitted_obj = splitted_obj.map{ object, clusters,  f ->
              def size = f.text.trim()
              tuple(object, clusters, size)
          }
      splitted_obj=splitted_obj.filter{it -> it[2].toInteger() >= 2*(params.minSize+1)} 
      splitted_obj=splitted_obj.map{t -> tuple(t[0],t[1])}
      splitted_obj.view()

      run_pipeline_iteration(splitted_obj)

      c1=run_pipeline_iteration.out.object.flatMap { list -> list } 
      c2=run_pipeline_iteration.out.clusters.flatMap { list -> list } 
      c3=run_pipeline_iteration.out.num_clusters.flatMap { list -> list } 
      ch=c1.merge(c2).merge(c3) 
      out_ch=ch.filter{it -> it[2].toInteger() > 1} 

    emit:
      out_ch=out_ch
}

// Run seventh iteration workflow
workflow run_iteration_recurse4{ 
    take:
     out_ch

     main:
      input_split=out_ch.map{t -> tuple(t[0],t[1],t[2].toInteger())}.flatMap{obj,cluster,num -> def nn=(1..num).toList()
                                             nn.collect{n -> tuple(obj,cluster,n)}}
      update_cluster_id(input_split)
      update_out=update_cluster_id.out.collect(flat: false).flatMap()
      input_split_updated = update_out.map{ object, clusters, cl, f ->
              def cl_updated = f.text.trim()
              tuple(object, clusters, cl, cl_updated)
          }
      split_object(input_split_updated,params.minSize) 
      splitted_obj = split_object.out.collect(flat: false).flatMap()
      splitted_obj.view()
      splitted_obj = splitted_obj.map{ object, clusters,  f ->
              def size = f.text.trim()
              tuple(object, clusters, size)
          }
      splitted_obj=splitted_obj.filter{it -> it[2].toInteger() >= 2*(params.minSize+1)} 
      splitted_obj=splitted_obj.map{t -> tuple(t[0],t[1])}
      splitted_obj.view()

      run_pipeline_iteration(splitted_obj)

      c1=run_pipeline_iteration.out.object.flatMap { list -> list } 
      c2=run_pipeline_iteration.out.clusters.flatMap { list -> list } 
      c3=run_pipeline_iteration.out.num_clusters.flatMap { list -> list } 
      ch=c1.merge(c2).merge(c3) 
      out_ch=ch.filter{it -> it[2].toInteger() > 1} 

    emit:
      out_ch=out_ch
}

// Run eigth iteration workflow
workflow run_iteration_recurse5{ 
    take:
     out_ch

     main:
      input_split=out_ch.map{t -> tuple(t[0],t[1],t[2].toInteger())}.flatMap{obj,cluster,num -> def nn=(1..num).toList()
                                             nn.collect{n -> tuple(obj,cluster,n)}}
      update_cluster_id(input_split)
      update_out=update_cluster_id.out.collect(flat: false).flatMap()
      input_split_updated = update_out.map{ object, clusters, cl, f ->
              def cl_updated = f.text.trim()
              tuple(object, clusters, cl, cl_updated)
          }
      split_object(input_split_updated,params.minSize) 
      splitted_obj = split_object.out.collect(flat: false).flatMap()
      splitted_obj.view()
      splitted_obj = splitted_obj.map{ object, clusters,  f ->
              def size = f.text.trim()
              tuple(object, clusters, size)
          }
      splitted_obj=splitted_obj.filter{it -> it[2].toInteger() >= 2*(params.minSize+1)} 
      splitted_obj=splitted_obj.map{t -> tuple(t[0],t[1])}
      splitted_obj.view()

      run_pipeline_iteration(splitted_obj)

      c1=run_pipeline_iteration.out.object.flatMap { list -> list } 
      c2=run_pipeline_iteration.out.clusters.flatMap { list -> list } 
      c3=run_pipeline_iteration.out.num_clusters.flatMap { list -> list } 
      ch=c1.merge(c2).merge(c3) 
      out_ch=ch.filter{it -> it[2].toInteger() > 1} 

    emit:
      out_ch=out_ch
}

workflow reassign_cells{
    take:
      out_ch
    
    main:
    reassign_discarded_cells(out_ch,params.out_path,params.tmp_path)
}

workflow {
    initialize(params.out_path, params.tmp_path)
    run_pipeline(params.input_object)
      out_ch = run_pipeline.out.object
        .combine(run_pipeline.out.clusters)
        .combine(run_pipeline.out.num_clusters)
      out_ch_filt=out_ch.filter{it -> it[2].toInteger() > 1}
    run_iteration(out_ch_filt)
    run_iteration_recurse(run_iteration.out)
    /*run_iteration_recurse.recurse(run_iteration.out).until{it -> it[0] > 1} recurse do not work with collect() Recursive workflows cannot use reduction operators such as collect, reduce, and toList, because these operators cause the recursion to hang indefinitely after the initial iteration. */
    run_iteration_recurse1(run_iteration_recurse.out)
    run_iteration_recurse2(run_iteration_recurse1.out)
    run_iteration_recurse3(run_iteration_recurse2.out)
    run_iteration_recurse4(run_iteration_recurse3.out)
    run_iteration_recurse5(run_iteration_recurse4.out)
      last_nonempty = run_iteration_recurse5.out.first().map{it -> it[2]}.ifEmpty { '1' }
    reassign_cells(last_nonempty)
}

