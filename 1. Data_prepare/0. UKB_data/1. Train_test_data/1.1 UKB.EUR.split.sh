library(data.table)

for (pop in c("EUR")) {
    # Load population data
    pop_id = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/", pop, "_id.tsv"))
    
    # Randomly shuffle indices
    set.seed(42)  # Ensure reproducibility
    total_n = nrow(pop_id)
    train_size = round(2/3 * total_n)  # Compute 2/3 size
    
    # Randomly sample indices
    train_indices = sample(seq_len(total_n), size = train_size)
    test_indices = setdiff(seq_len(total_n), train_indices)  # Remaining 1/3

    # Split into training and testing
    pop_train_id = pop_id[train_indices, ]
    pop_test_id = pop_id[test_indices, ]

    # Save outputs
    write.table(pop_train_id, paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/", pop, "_train_2over3_doubleid.tsv"), 
                row.names=F, col.names=F, quote = F, sep = "\t")
    
    write.table(pop_test_id, paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/", pop, "_test_1over3_doubleid.tsv"),
                row.names=F, col.names=F, quote = F, sep = "\t")
}

## full EUR: /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_id.tsv
## train EUR: /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/ancestry_info/EUR_train_2over3_doubleid.tsv