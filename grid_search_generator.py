import os
import json
"""Script that generates files for grid search configs"""

gcn_configs = os.listdir("configs")
#batch_sizes = [20, 35, 50, 65, 80, 95]
batch_sizes = [500, 1000]
neighbors = [10, 50, 100, 250]
#dropouts = [0.0, 0.2, 0.4, 0.6, 0.8]
dropouts = [0.0]

counter = 1
for config in gcn_configs:
    for dropout in dropouts:
        for batch in batch_sizes:
            for n in neighbors:
                if n >= batch: continue
                d = {"config":config, "dropout":dropout, "batch_size": batch, "neighbors":n}

                with open(("grid_search_files_103022/" + str(counter)+ ".txt"), 'w') as output:
                    output.write(json.dumps(d))
                    
                counter +=1

print(counter)
