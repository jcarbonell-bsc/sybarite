
scFEA_wrapper <- function(norm_counts, scfea_path){
  
  tmp <- tempfile(pattern = "tmp_scfea__", tmpdir = ".")
  dir.create(tmp)
  
  tmp_file <- paste0(tmp, "/norm_counts.csv")
  write.table(norm_counts, file=tmp_file, row.names=T, col.names = T, sep=",", quote = T)
  
  command <- paste0("python ", scfea_path, "/src/scFEA.py --data_dir ", scfea_path, "/data --input_dir ", tmp, " --res_dir ", tmp, " --output_flux_file ", tmp, "/flux.csv --output_balance_file ", tmp, "/balance.csv --test_file ", basename(tmp_file))
  print(command)
  
  system(command)
  
  flux <- t(read.table(paste0(tmp, "/flux.csv"), sep=",", header=T, row.names = 1))
  balance <- t(read.table(paste0(tmp, "/balance.csv"), sep=",", header=T, row.names = 1))
  
  return(list(
    flux=flux,
    balance=balance
  ))
  
}