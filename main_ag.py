# Import sulla directory delle librerie
import sys
sys.path.append('libs')
from GeneticAlgorithm import GeneticAlgorithm
# Librerie di utilita
import csv
import math
import matplotlib.pyplot as plt
import numpy



# Comparazione della fitness
def fitness_cmp_max(val_1, val_2):
    return 1 if val_1 > val_2 else 0


def fitness_cmp_min(val_1, val_2):
    return 1 if val_1 < val_2 else 0


def compute_scores(mirna_data, up_c):
    scores = [[0 for i in range(2)] for j in range(len(mirna_data))]
    for k in range(len(mirna_data)):
        scores[k][0] = sum(mirna_data[k][0:up_c])
        scores[k][1] = sum(mirna_data[k][up_c:])
    return scores


def compute_scores2(mirna_data, up_c):
    scores_pos = [[] for j in range(len(mirna_data))]
    scores_neg = [[] for j in range(len(mirna_data))]
    for i in range(len(mirna_data)):
        for k in range(up_c):
            if mirna_data[i][k] == 1:
                scores_pos[i].append(k)
        for l in range(up_c + 1, len(mirna_data[i])):
            if mirna_data[i][l] == 1:
                scores_neg[i].append(l)
    return scores_pos, scores_neg


# Funzione/i Obiettivo
def fitness1(solution):
    return sum(solution)


def fitness2(solution):
    pos = 0.0
    neg = 0.0
    utl = 0.0
    gamma = 2
    eps = 3
    psi = 3
    for i in range(len(solution)):
        if solution[i] == 1:
            pos += scores[i][0]
            neg += scores[i][1]
            utl += float(scores[i][0]) / (float(scores[i][1]) + 1)
    mirnas = sum(solution)
    final_score = (float(math.pow(pos, gamma))) / (float(math.pow(mirnas, eps)) * float(math.pow(neg, psi)))
    return final_score


def fitness3(solution):
    pos = set() # np_regulated impattate
    neg = set() # normali impattate
    for i in range(len(solution)):
        if solution[i] == 1:
            for j in range(len(scores_pos[i])):
                pos.add(scores_pos[i][j])
            for k in range(len(scores_neg[i])):
                neg.add(scores_neg[i][k])
    # Il numero di elementi dei set mi dici quanto sono state impattate le proteine
    # dunque quanto sono distribuiti gli 1 li in mezzo.
    B1 = 10.0
    B2 = (1.0 / float(len(mirna_data[0]) - up_reg_cnt)) * up_reg_cnt
    final_score = ((B1 * math.fabs(len(pos) - up_reg_cnt)) + (B2 * (len(neg)))) * math.pow(sum(solution), 1)
    return final_score


# Lettura dei dati
def data_read():
    protein_classes = {}
    up_reg_cnt = 0
    # Lettura delle classi
    with open('dataset/proteins_classes_BRCA_ER-.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(csvreader) # Skippa la prima riga che rappresenta header
        for row in csvreader:
            # Aggiungo identificativo di classe
            # 1 = up_regulated, 0 = Altrimenti
            if int(row[1]) == 1:
                protein_classes[row[0]] = 'U'
                up_reg_cnt += 1
            else:
                protein_classes[row[0]] = 'N'
    # Lettura della matrice del secondo file
    with open('dataset/prot_mirna_matrix_L0_BRCA_ER-.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        dataset = []
        mirna_id = next(csvreader)[1:]
        for row in csvreader:
            row.insert(0, protein_classes[row[0]]) # Aggiungo etichetta di classe della proteina
            dataset.append(row)
    # Riordino le righe del dataset in base all'ordinamento delle classi di proteine
    # Nella prima parte ho tutte le up_regolated, nella seconda tutte le normali
    dataset = sorted(dataset, reverse=True, key=lambda row: row[0])
    mirna_data = []
    # salvo i dati dei singoli miRNA
    for i in range(2, len(dataset[0])):
        mirna_data.append([int(dataset[j][i]) for j in range(0, len(dataset))])
    return dataset, mirna_id, mirna_data, up_reg_cnt


def save_solution(solution, filepath):
    # results/metodo_mio/solutions.csv
    # results/metodo_prof/solutions.csv
    file = open(filepath, "a")
    string = ""
    for bit in solution:
        string += str(bit)
    file.write(string + "\n")
    file.close()


# Calcola la fequenza di utilizzo dei mirna nelle soluzioni trovate
def compute_mirna_frequency(filepath):
    # results/metodo_mio/solutions.csv
    # results/metodo_mio/solutions.csv
    with open(filepath, "r") as file:
        lines = file.readlines()
        sols = []
        for line in lines:
            sol = []
            for i in range(len(line) - 1):
                sol.append(int(line[i]))
            sols.append(sol)
        # calcolo delle frequenze
        frq = [sum(sols[i][j] for i in range(len(sols))) for j in range(len(sols[0]))]
    return frq


# Prendo i primi n migliori / peggiori
def get_bst_wst(n, frq):
    std_frq_idx = numpy.argsort(frq)
    wst = std_frq_idx[0:n]
    bst = std_frq_idx[-n-1:-1]
    return bst, wst


# Carica i dati dei miei run
def load_run_data(path):
    # results/metodo_mio/runs/
    # results/metodo_prof/runs/
    all_dataset = []
    for i in range(100):
        with open(path + "run_"+str(i) + ".csv") as file:
            csvreader = csv.reader(file, delimiter=';', quotechar='|')
            dataset = []
            for row in csvreader:
                dataset.append([int(v) for v in row])
            all_dataset.append(dataset)
    # Adesso dobbiamo normalizzare i datasets
    # vediamo chi e' dataset piu' lunggo
    longest = 0
    for i in range(len(all_dataset)):
        if len(all_dataset[i]) > longest:
            longest = len(all_dataset[i])
    # una volta capito l'indice del dataset piu' lungo
    # aggiungiamo righe a quelli piu' corti affiche
    # tutti abbiano la stessa dimensione
    for i in range(len(all_dataset)):
        # se e' piu' corto del piu' grande
        # allora aggiungo righe
        if len(all_dataset[i]) < longest:
            lst_elm = all_dataset[i][-1]
            for j in range(len(all_dataset[i]), longest):
                all_dataset[i].append(lst_elm)
    # calcolo del dataset unico normalizzato
    all_ds_nrm = [[0 for i in range(3)] for j in range(longest)]
    # adesso posso tornare il dataset normalizzato
    for x in range(len(all_dataset[0])):
        # riga x-sima presa da tutti i dataset
        rows = [all_dataset[i][x] for i in range(len(all_dataset))]
        all_ds_nrm[x][0] = sum(rows[i][0] for i in range(len(rows)))
        all_ds_nrm[x][1] = sum(rows[i][1] for i in range(len(rows)))
        all_ds_nrm[x][2] = sum(rows[i][2] for i in range(len(rows)))
    all_ds_nrm = [[float(all_ds_nrm[i][j])/float(len(all_dataset)) for j in range(len(all_ds_nrm[0]))] for i in range(len(all_ds_nrm))]
    return all_ds_nrm



def plot_mirna(frq, bst, wst):
    lab = [mirna_id[bst[i]] for i in range(len(bst))] + [mirna_id[wst[i]] for i in range(len(wst))]
    wgt = [float(frq[bst[i]])/float(sum(frq)) for i in range(len(bst))] + [float(frq[wst[i]])/float(sum(frq)) for i in range(len(wst))]
    all = list(bst) + list(wst)
    fig, axs = plt.subplots(4, 1, figsize=(16, 9), shared=False)
    # istogramma delle frequenze
    axs[0].set_ylabel("frequency")
    axs[0].bar(lab, wgt)
    # grafico delle utilita
    pos = [float(scores[i][0]) for i in all]
    axs[1].set_ylabel("up_proteins")
    axs[1].plot(lab, pos)
    neg = [float(scores[i][1]) for i in all]
    axs[2].set_ylabel("nr_proteins")
    axs[2].plot(lab, neg)
    utl = [float(scores[i][0])/float(scores[i][1]) for i in all]
    axs[3].set_ylabel("up_proteins/nr_proteins")
    axs[3].plot(lab, utl)
    # Rotazione della label
    for label in axs[0].get_xmajorticklabels() + axs[1].get_xmajorticklabels() + axs[2].get_xmajorticklabels() + axs[3].get_xmajorticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment("right")
    plt.show()



# Plotta le performance dell'algoritmo genetico
# dando in ingresso il dataset normalizzato
def plot_avg_ag_perf(nrm_dts):
    fig = plt.figure()
    times = range(len(nrm_dts))
    np_covered = [nrm_dts[i][0] for i in range(len(nrm_dts))]
    nr_covered = [nrm_dts[i][1] for i in range(len(nrm_dts))]
    mirnas = [nrm_dts[i][2] for i in range(len(nrm_dts))]

    plt.subplot(2, 2, 1)
    plt.xlabel('generations')
    plt.ylabel('np_covered')
    plt.plot(times, np_covered)
    plt.grid(True)

    plt.subplot(2, 2, 2)
    plt.xlabel('generations')
    plt.ylabel('nr_covered')
    plt.plot(times, nr_covered)
    plt.grid(True)

    plt.subplot(2, 2, 3)
    plt.xlabel('generations')
    plt.ylabel('subsets')
    plt.plot(times, mirnas)
    plt.grid(True)

    plt.show()


# Main Program
if __name__ == "__main__":

    # lettura e caricamento dati
    dataset, mirna_id, mirna_data, up_reg_cnt = data_read()
    scores = compute_scores(mirna_data, up_reg_cnt)
    scores_pos, scores_neg = compute_scores2(mirna_data, up_reg_cnt)

    # effettuo analisi sui mirna, e sulla loro frequenza
    # di prelevamento in base alla performance
    '''
    frq = compute_mirna_frequency("results/metodo_mio/solutions.csv")
    bst, wst = get_bst_wst(10, frq)
    plot_mirna(frq, bst, wst)
    '''

    # Normalizzo i dataset e ne ottendo un'unico normalizzato
    '''
    all_dataset_norm = load_run_data("results/metodo_mio/runs/")
    plot_avg_ag_perf(all_dataset_norm)
    '''

    # esecuzione dell'algoritmo genetico, su un numero elevato di run
    for i in range(100):
        sol_len = len(mirna_data)
        GA = GeneticAlgorithm(fitness3, fitness_cmp_min, 0, 1, 100, sol_len)
        GA.set_params(mirna_data, up_reg_cnt)
        sol = GA.run(gen=20000, pc=0.8, pm=2.0/sol_len, debug=True)
        GA.save_data(i, "results/metodo_mio/runs/")
        save_solution(sol,"results/metodo_mio/solutions.csv")