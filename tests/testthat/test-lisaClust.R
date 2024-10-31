test_that(
  "the output is equal to original result", 
  {
    load(system.file("testdata/original_result.rda", package = "lisaClust"))
    set.seed(51773)
    
    kerenSPE <- SpatialDatasets::spe_Keren_2018()
    kerenSPE <- kerenSPE[,kerenSPE$imageID %in% c("5", "6")]
    kerenSPE <- suppressWarnings(lisaClust(kerenSPE, k = 5)) |>
      colData()
    
    expect_equal(
      suppressWarnings(kerenSPE),
      original_result
    )
  }
)
