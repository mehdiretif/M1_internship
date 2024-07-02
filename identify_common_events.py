#!/usr/bin/env python

def reading_file(path):
    file=open(path, 'r')
    text=file.readlines()
    file.close()
    return text

def parsing_se_rmats(path):
    dico_se={}
    file=reading_file(path)
    for i, line in enumerate(file):
        if i ==0:
            continue
        elements = line.strip().split("\t")
        gene = elements[2].strip('"')
        chr = elements[3]
        A = int(elements[5])
        B = int(elements[6])
        event_ID = elements[0]
        vast_type_line=f"{chr}:{A}-{B}"

        if gene not in dico_se:
            dico_se[gene] = {}

        if event_ID not in dico_se[gene]:
            dico_se[gene][event_ID] = []
        dico_se[gene][event_ID].append(vast_type_line)

    return dico_se


def vast_to_dico(path_vast):
    vast_file=reading_file(path_vast)
    dico_ex={}
    for i, line in enumerate(vast_file):
        if i ==0:
            continue
        
        elements = line.strip().split("\t")
        gene=elements[0]
        interest=elements[6].split(":")
        chr = interest[0]
        coord = interest[1].split("-")
        A = int(coord[0])
        B = int(coord[1])
        event_ID = elements[1]
       
        if gene not in dico_ex:
            dico_ex[gene] = {}
        if event_ID not in dico_ex[gene]:
            dico_ex[gene][event_ID] = []

        for deltaA in range(-5,6):
            for deltaB in range(-5,6):
                new_A = A + deltaA
                new_B = B + deltaB
                new_coord=f"{chr}:{new_A}-{new_B}"
                dico_ex[gene][event_ID].append(new_coord)
                
    return len(dico_ex['GLI4']['HsaEX6082696'])

def common_event_research(path_vast, path_rmats):
    dico_se=parsing_se_rmats(path_rmats)
    dico_ex=vast_to_dico(path_vast)
    common_events={}
    for gene in dico_ex.keys():  
        if gene in dico_se.keys():
            for rmats_eventID, rmats_coord in dico_se[gene].items():
                #print(rmats_coord)
                for vast_eventID, all_vast_coords in dico_ex[gene].items():
                    for vast_coord in all_vast_coords:
                        #print(vast_coord)
                        if rmats_coord[0]==vast_coord:
                            if gene not in common_events:
                                common_events[gene] = []
                            common_events[gene].append([rmats_eventID, vast_eventID, rmats_coord[0]])  
    return(common_events)      
                            
def create_common_rmats_vast_file(path_vast, path_rmats, output):
    common_dico =  common_event_research(path_vast, path_rmats)       
    new_file = open(output, 'w')
    new_file.write("Gene" + "\t" + "rMats_event_ID" + "\t" + "Vast_ID" + "\t" + "Exon_coordinates" + "\n")
    for gene, content in common_dico.items():
        for event_info in content:
            new_file.write(gene +"\t" + event_info[0] + "\t" + event_info[1] + "\t" + event_info[2] + "\n")
    new_file.close()
                

            

#print(parsing_se_rmats("rmats/SE_significant.MATS.JC.txt"))
print(vast_to_dico("../vast/results/nf_results/EX_significant.txt"))
#print(common_event_research("vast/results/nf_results/EX_significant.txt", "rmats/SE_significant.MATS.JC.txt"))

#create_common_rmats_vast_file("vast/results/nf_results/EX_raw_inclusion_file.DIFF.txt", "rmats/SE_significant.MATS.JC.txt","common_events_based_non_significant_vast.txt")

#create_common_rmats_vast_file("vast/results/nf_results/EX_significant.txt", "rmats/SE_significant.MATS.JC.txt","common_events.txt")


