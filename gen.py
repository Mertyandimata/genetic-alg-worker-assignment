import random
import copy
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box

# Sabit random seed: Her çalıştırmada aynı rastgele veriler üretilsin
random.seed(42)
console = Console()

# =============================================================================
# 1. WORKCENTER DATABASE
# =============================================================================
production_lines = [f"Line{i}" for i in range(1, 5)]
workcenters = []
wc_id = 1
for line in production_lines:
    num_wcs = random.randint(4, 7)
    for _ in range(num_wcs):
        wc = {
            "ProductionLineName": line,
            "WorkcenterName": f"WC{wc_id}",
            "Capacity": random.randint(4, 8),
            "MinTrainingLevel": random.choice([1, 2]),
            "Q_Limitation1": random.choice([0,5, 10]),
            "HighLiftRequired": random.choice([0, 1])
        }
        workcenters.append(wc)
        wc_id += 1
wc_map = {wc["WorkcenterName"]: wc for wc in workcenters}

# =============================================================================
# 2. WORKER DATABASE
#
# Her çalışanın:
#  - WorkerID, Name
#  - Limitation: %90 = 0 (limitsiz), %5 = 5, %5 = 10
#  - HighLiftCapability: %10 olasılıkla 1, %90 = 0
#  - CandidateWCs: Çalışanın atanabileceği workcenter isimlerinden (global wc_map'den)
#                  rastgele 5-8 tanesi seçilir ve her biri için ayrı training seviyesi (1,2,3) atanır.
# =============================================================================
NUM_WORKERS = 60
worker_db = []
all_wc_names = list(wc_map.keys())
for i in range(1, NUM_WORKERS + 1):
    p = random.random()
    if p < 0.90:
        lim = 0
    elif p < 0.95:
        lim = 5
    else:
        lim = 10
    p2 = random.random()
    high_lift = 1 if p2 < 0.10 else 0
    num_candidates = random.randint(5, 8)
    candidates = random.sample(all_wc_names, num_candidates)
    candidate_dict = {wc: random.choice([ 1,2,3]) for wc in candidates}
    worker = {
        "WorkerID": i,
        "Name": f"Worker{i}",
        "Limitation": lim,
        "HighLiftCapability": high_lift,
        "CandidateWCs": candidate_dict
    }
    worker_db.append(worker)

# =============================================================================
# 3. DATABASE VISUALIZATION FUNCTIONS
# =============================================================================
def print_worker_db():
    table = Table(title="Worker Database", box=box.MINIMAL_DOUBLE_HEAD)
    table.add_column("WID", justify="right")
    table.add_column("Name")
    table.add_column("Limit", justify="center")
    table.add_column("HighLift", justify="center")
    table.add_column("Candidate WCs (TrainLv)", justify="left")
    for w in worker_db:
        cand_str = ", ".join(f"{wc}({tr})" for wc, tr in w["CandidateWCs"].items())
        table.add_row(str(w["WorkerID"]), w["Name"], str(w["Limitation"]),
                      str(w["HighLiftCapability"]), cand_str)
    console.print(table)

def print_workcenter_db():
    table = Table(title="Workcenter Database", box=box.MINIMAL_DOUBLE_HEAD)
    table.add_column("Line")
    table.add_column("WC Name", justify="center")
    table.add_column("Cap", justify="right")
    table.add_column("MinTrLv", justify="center")
    table.add_column("Q_Lim(kg)", justify="center")
    table.add_column("HighLiftReq", justify="center")
    for wc in workcenters:
        table.add_row(wc["ProductionLineName"], wc["WorkcenterName"],
                      str(wc["Capacity"]), str(wc["MinTrainingLevel"]),
                      str(wc["Q_Limitation1"]), str(wc["HighLiftRequired"]))
    console.print(table)

print_worker_db()
print_workcenter_db()

# =============================================================================
# 4. TRY_REASSIGN FUNCTION
#
# Çalışanın mevcut ataması uygun değilse, çalışanın CandidateWCs listesindeki WC'ler arasında
# (mevcut atama, kapasite durumu göz önünde bulundurularak) uygun bir alternatif arar.
# =============================================================================
def try_reassign(worker_idx, current_assignment):
    w = worker_db[worker_idx]
    candidates = list(w["CandidateWCs"].keys())
    random.shuffle(candidates)
    for cand in candidates:
        wc_candidate = wc_map[cand]
        # Eğitim seviyesi kontrolü eklendi
        if w["CandidateWCs"][cand] < wc_candidate["MinTrainingLevel"]:
            continue
        if w["Limitation"] != 0 and w["Limitation"] < wc_candidate["Q_Limitation1"]:
            continue
        if wc_candidate["HighLiftRequired"] == 1 and w["HighLiftCapability"] != 1:
            continue
        count = sum(1 for j, ass in enumerate(current_assignment) 
                    if j != worker_idx and ass == cand)
        if count < wc_candidate["Capacity"]:
            return cand, None
    return None, "No feasible alternative found"

# =============================================================================
# 5. IMPROVED REPAIR FUNCTION (GÜNCELLENDİ)
# =============================================================================
def improved_repair(individual):
    max_iter = 3
    current_assignment = copy.deepcopy(individual)
    unassigned_reasons = {}
    for _ in range(max_iter):
        changed = False
        for i, wc_name in enumerate(current_assignment):
            w = worker_db[i]
            if wc_name is not None and wc_name not in w["CandidateWCs"]:
                new_wc, reason = try_reassign(i, current_assignment)
                if new_wc is not None:
                    current_assignment[i] = new_wc
                    changed = True
                else:
                    current_assignment[i] = None
                    unassigned_reasons[w["WorkerID"]] = "Assigned WC not in candidate list and no alternative"
                continue
            if wc_name is None:
                new_wc, reason = try_reassign(i, current_assignment)
                if new_wc is not None:
                    current_assignment[i] = new_wc
                    changed = True
                else:
                    unassigned_reasons[w["WorkerID"]] = reason
                continue
            if wc_name not in wc_map:
                new_wc, reason = try_reassign(i, current_assignment)
                if new_wc is not None:
                    current_assignment[i] = new_wc
                    changed = True
                else:
                    current_assignment[i] = None
                    unassigned_reasons[w["WorkerID"]] = "Invalid WC and no alternative"
                continue
            wc = wc_map[wc_name]
            if w["Limitation"] != 0 and w["Limitation"] < wc["Q_Limitation1"]:
                new_wc, reason = try_reassign(i, current_assignment)
                if new_wc is not None:
                    current_assignment[i] = new_wc
                    changed = True
                else:
                    current_assignment[i] = None
                    unassigned_reasons[w["WorkerID"]] = f"Worker limit {w['Limitation']} < required {wc['Q_Limitation1']} and no alternative"
                continue
            # Eğitim seviyesi direkt olarak CandidateWCs'den alınıyor
            candidate_training = w["CandidateWCs"][wc_name]  # .get() KALDIRILDI
            if candidate_training < wc["MinTrainingLevel"]:
                new_wc, reason = try_reassign(i, current_assignment)
                if new_wc is not None:
                    current_assignment[i] = new_wc
                    changed = True
                else:
                    current_assignment[i] = None
                    unassigned_reasons[w["WorkerID"]] = f"Worker training {candidate_training} < required {wc['MinTrainingLevel']} and no alternative"
                continue
            if wc["HighLiftRequired"] == 1 and w["HighLiftCapability"] != 1:
                new_wc, reason = try_reassign(i, current_assignment)
                if new_wc is not None:
                    current_assignment[i] = new_wc
                    changed = True
                else:
                    current_assignment[i] = None
                    unassigned_reasons[w["WorkerID"]] = f"WC {wc_name} requires high lift but worker lacks capability"
        # Kapasite kontrolü aynı kaldı
        wc_assign = {}
        for i, ass in enumerate(current_assignment):
            if ass is not None:
                wc_assign.setdefault(ass, []).append(i)
        for wcn, idx_list in wc_assign.items():
            cap = wc_map[wcn]["Capacity"]
            if len(idx_list) > cap:
                idx_list.sort(key=lambda ix: worker_db[ix]["CandidateWCs"][wcn])  # .get() KALDIRILDI
                excess = len(idx_list) - cap
                for ix in idx_list[:excess]:
                    new_wc, reason = try_reassign(ix, current_assignment)
                    if new_wc is not None:
                        current_assignment[ix] = new_wc
                        changed = True
                    else:
                        current_assignment[ix] = None
                        unassigned_reasons[worker_db[ix]["WorkerID"]] = f"Over capacity in {wcn} and no alternative"
        if not changed:
            break
    return current_assignment, unassigned_reasons

# =============================================================================
# 6A. FORCE TRAINING FUNCTIONALITY
#
# Eğer çalışanın aday training değeri, atandığı WC'nin MinTrainingLevel'ından düşükse,
# forced training cost (-200) uygulanır.
# =============================================================================
def forced_training_cost(worker_idx, current_assignment):
    w = worker_db[worker_idx]
    if current_assignment[worker_idx] is None:
        return 0
    wc = wc_map[current_assignment[worker_idx]]
    if w["CandidateWCs"].get(current_assignment[worker_idx], 0) < wc["MinTrainingLevel"]:
        return -200
    return 0

# =============================================================================
# 7. FITNESS FUNCTION
#
# - Atanan her çalışan için:
#     Base: +100
#     Training:
#       * Eğer aday training < WC.MinTrainingLevel: effective as equal (+20) but forced cost (-200) uygulanır.
#       * Eğer eşitse: +20
#       * Eğer üstündeyse: +20 + (diff * 30)
# - Eğer worker.Limitation == 0: +10 bonus
# - Full kapasite WC: +50 bonus
# - Unassigned her çalışan için: -150 ceza
# - Tüm çalışan atandıysa: +1000 bonus
# =============================================================================
def evaluate_individual(individual):
    repaired, reasons = improved_repair(individual)
    score = 0
    for wcn in wc_map:
        cap = wc_map[wcn]["Capacity"]
        assigned = sum(1 for x in repaired if x == wcn)
        if assigned == cap and cap > 0:
            score += 50
    unassigned_count = 0
    for i, wcn in enumerate(repaired):
        if wcn is None:
            score -= 150
            unassigned_count += 1
        else:
            w = worker_db[i]
            wc = wc_map[wcn]
            base = 100
            candidate_training = w["CandidateWCs"][wcn]  # .get() KALDIRILDI
            if candidate_training < wc["MinTrainingLevel"]:
                base += 20
                base += forced_training_cost(i, repaired)
            else:
                diff = candidate_training - wc["MinTrainingLevel"]
                base += 20 + diff * 30
            if w["Limitation"] == 0:
                base += 10
            score += base
    if unassigned_count == 0:
        score += 1000
    return score, repaired, reasons

# =============================================================================
# 8. SMART INITIAL POPULATION
#
# Her çalışan için, yalnızca kendi CandidateWCs listesinden seçim yapılır.
# =============================================================================
def generate_initial_population_smart(pop_size):
    population = []
    for _ in range(pop_size):
        indiv = []
        for w in worker_db:
            candidates = list(w["CandidateWCs"].keys())
            if candidates:
                indiv.append(random.choice(candidates))
            else:
                indiv.append(None)
        population.append(indiv)
    return population

# =============================================================================
# 9. GA OPERATORS: Selection, 2-Point Crossover, Mutation, Random Immigrants
# =============================================================================
def selection(population, fitnesses, k=2):
    parents = []
    t_size = 4
    pop_size = len(population)
    for _ in range(k):
        cands = random.sample(range(pop_size), t_size)
        winner = max(cands, key=lambda c: fitnesses[c])
        parents.append(copy.deepcopy(population[winner]))
    return parents

def crossover_2point(p1, p2):
    size = len(p1)
    c1, c2 = copy.deepcopy(p1), copy.deepcopy(p2)
    cut1, cut2 = sorted(random.sample(range(size), 2))
    
    # Çaprazlama sonrası geçersiz atamaları düzelt
    for i in range(cut1, cut2):
        w = worker_db[i]
        if c1[i] not in w["CandidateWCs"]:
            c1[i] = random.choice(list(w["CandidateWCs"].keys())) if w["CandidateWCs"] else None
        if c2[i] not in w["CandidateWCs"]:
            c2[i] = random.choice(list(w["CandidateWCs"].keys())) if w["CandidateWCs"] else None
            
    return c1, c2

def mutation(indiv, m_rate=0.1):
    for i in range(len(indiv)):
        if random.random() < m_rate:
            w = worker_db[i]
            candidates = list(w["CandidateWCs"].keys())
            if candidates:
                indiv[i] = random.choice(candidates)  # SADECE CandidateWCs'den seç
            else:
                indiv[i] = None
    return indiv

def generate_random_individual():
    indiv = []
    for w in worker_db:
        candidates = list(w["CandidateWCs"].keys())
        if candidates:
            indiv.append(random.choice(candidates))  # SADECE CandidateWCs'den seç
        else:
            indiv.append(None)
    return indiv

# =============================================================================
# 10. GENETIC ALGORITHM (ADVANCED)
#
# - Popülasyonun %5’i random immigrants ile yenilenir.
# - Erken durdurma:
#     a) En iyi bireyde hiç None yoksa (tüm çalışan atandıysa),
#     b) veya stagnation limitine ulaşırsa.
# =============================================================================
def genetic_algorithm(pop_size=50, generations=1000, mutation_rate=0.1,
                      elitism=2, stagnation_limit=200, random_immigrants_rate=0.05):
    population = generate_initial_population_smart(pop_size)
    best_fit = float('-inf')
    best_indiv = None
    best_rep = None
    best_unassigned = {}
    stagnation = 0
    for gen in range(generations):
        fitnesses = []
        rep_list = []
        unass_list = []
        for indiv in population:
            sc, rep, unass = evaluate_individual(indiv)
            fitnesses.append(sc)
            rep_list.append(rep)
            unass_list.append(unass)
        max_fit = max(fitnesses)
        if max_fit > best_fit:
            best_fit = max_fit
            idx = fitnesses.index(max_fit)
            best_indiv = population[idx]
            best_rep = rep_list[idx]
            best_unassigned = unass_list[idx]
            stagnation = 0
        else:
            stagnation += 1
        none_count = sum(1 for x in best_rep if x is None)
        if none_count == 0:
            console.log(f"All workers assigned at Gen {gen+1}. Early stopping.")
            break
        if stagnation >= stagnation_limit:
            console.log(f"Early stopping at Gen {gen+1} (no improvement for {stagnation_limit} gens).")
            break
        new_population = []
        sorted_idx = sorted(range(len(population)), key=lambda i: fitnesses[i], reverse=True)
        for s in sorted_idx[:elitism]:
            new_population.append(copy.deepcopy(population[s]))
        while len(new_population) < pop_size:
            p = selection(population, fitnesses, 2)
            c1, c2 = crossover_2point(p[0], p[1])
            c1 = mutation(c1, mutation_rate)
            c2 = mutation(c2, mutation_rate)
            new_population.append(c1)
            if len(new_population) < pop_size:
                new_population.append(c2)
        num_immigrants = int(random_immigrants_rate * pop_size)
        for _ in range(num_immigrants):
            idx = random.randint(elitism, pop_size - 1)
            new_population[idx] = generate_random_individual()
        population = new_population
        console.log(f"Gen {gen+1} | BestFit = {best_fit} | Unassigned = {none_count}")
    return best_rep, best_fit, best_unassigned

# =============================================================================
# 11. LOCAL SWAP IMPROVEMENT
#
# GA sonrasında, çiftler arası swap yaparak fitness artışı aranır.
# =============================================================================
def local_swap_improvement(assignment):
    best_assignment = copy.deepcopy(assignment)
    best_score, _, _ = evaluate_individual(best_assignment)
    improved = True
    while improved:
        improved = False
        for i in range(len(best_assignment)):
            for j in range(i+1, len(best_assignment)):
                # Swap işlemi sadece aynı CandidateWCs'deki WC'ler arasında
                if (best_assignment[i] in worker_db[j]["CandidateWCs"] and 
                    best_assignment[j] in worker_db[i]["CandidateWCs"]):
                    candidate = copy.deepcopy(best_assignment)
                    candidate[i], candidate[j] = candidate[j], candidate[i]
                    candidate_score, _, _ = evaluate_individual(candidate)
                    if candidate_score > best_score:
                        best_assignment = candidate
                        best_score = candidate_score
                        improved = True
                        break
            if improved:
                break
    return best_assignment, best_score

# =============================================================================
# 12. DETAILED UNASSIGNED REPORT
#
# Her unassigned (None) çalışan için, çalışanın özelliklerini, atanamama nedenini ve
# aday listesindeki alternatif WC'leri raporlar.
# =============================================================================
def find_feasible_alternatives(worker):
    # Çalışanın CandidateWCs listesindeki tüm workcenter'ları döndürür.
    return list(worker["CandidateWCs"].keys())

def print_detailed_unassigned_report(unassigned_info):
    table = Table(title="Detailed Unassigned Report", box=box.MINIMAL_DOUBLE_HEAD)
    table.add_column("WorkerID", justify="right")
    table.add_column("Name")
    table.add_column("TrainLv", justify="center")
    table.add_column("Limit", justify="center")
    table.add_column("HighLift", justify="center")
    table.add_column("Reason")
    table.add_column("Feasible Alternatives")
    for wid, reason in unassigned_info.items():
        worker = next((w for w in worker_db if w["WorkerID"] == wid), None)
        if worker is None:
            continue
        alts = find_feasible_alternatives(worker)
        alt_str = ", ".join(alts) if alts else "None"
        table.add_row(str(worker["WorkerID"]),
                      worker["Name"],
                      "-",  # Her WC için ayrı training seviyesi olduğundan genel bir değer verilemez.
                      str(worker["Limitation"]),
                      str(worker["HighLiftCapability"]),
                      reason,
                      alt_str)
    console.print(table)

# =============================================================================
# 13. GENERAL STATISTICS CALCULATION
# =============================================================================
def calculate_statistics(assignment):
    stats = {}
    total_workers = len(worker_db)
    assigned_workers = [w for i, w in enumerate(worker_db) if assignment[i] is not None]
    unassigned_workers = [w for i, w in enumerate(worker_db) if assignment[i] is None]
    stats["Total Workers"] = total_workers
    stats["Assigned Workers"] = len(assigned_workers)
    stats["Unassigned Workers"] = len(unassigned_workers)
    wc_stats = {}
    for wcn in wc_map:
        count = sum(1 for i, a in enumerate(assignment) if a == wcn)
        wc_stats[wcn] = {"Capacity": wc_map[wcn]["Capacity"], "Assigned": count, "Full": count == wc_map[wcn]["Capacity"]}
    stats["Workcenter Stats"] = wc_stats
    forced_training_count = 0
    for i, a in enumerate(assignment):
        if a is not None:
            w = worker_db[i]
            effective_train = w["CandidateWCs"].get(a, 0)
            wc = wc_map[a]
            if effective_train < wc["MinTrainingLevel"]:
                forced_training_count += 1
    stats["Forced Training Count"] = forced_training_count
    return stats

def print_statistics(stats):
    lines = []
    lines.append(f"Total Workers: {stats['Total Workers']}")
    lines.append(f"Assigned Workers: {stats['Assigned Workers']}")
    lines.append(f"Unassigned Workers: {stats['Unassigned Workers']}")
    lines.append(f"Forced Training Count: {stats['Forced Training Count']}")
    lines.append("Workcenter Utilization:")
    for wcn, data in stats["Workcenter Stats"].items():
        status = "Full" if data["Full"] else "Not Full"
        lines.append(f"  - {wcn}: {data['Assigned']}/{data['Capacity']} ({status})")
    console.print(Panel("\n".join(lines), title="General Statistics", expand=False))

# =============================================================================
# 14. GA MAIN
# =============================================================================
best_assignment, best_score, unassigned_info = genetic_algorithm()

# GA sonrasında local swap improvement
final_assignment, final_score = local_swap_improvement(best_assignment)
stats = calculate_statistics(final_assignment)

# =============================================================================
# 15. FINAL RESULT REPORTING
# =============================================================================
def print_final_result(assignment, score, unassigned_info):
    table = Table(title=f"Best Assignment (Score = {score})", box=box.MINIMAL_DOUBLE_HEAD)
    table.add_column("WorkerID", justify="right")
    table.add_column("Name")
    table.add_column("TrainLv", justify="center")
    table.add_column("Limit", justify="center")
    table.add_column("WC", justify="center")
    table.add_column("Line", justify="center")
    table.add_column("Status", justify="center")
    for i, wcn in enumerate(assignment):
        w = worker_db[i]
        if wcn is None:
            table.add_row(str(w["WorkerID"]), w["Name"],
                          "-", str(w["Limitation"]),
                          "-", "-",
                          "[red]UNASSIGNED[/red]")
        else:
            wc = wc_map[wcn]
            effective_train = w["CandidateWCs"][wcn]  # .get() KALDIRILDI
            status = "[green]OK[/green]" if effective_train >= wc["MinTrainingLevel"] else "[yellow]LOW TRAINING[/yellow]"
            table.add_row(str(w["WorkerID"]), w["Name"],
                          str(effective_train), str(w["Limitation"]),
                          wcn, wc["ProductionLineName"], status)
    console.print(table)
    console.print(Panel(f"Final Score: {score}", title="Final Fitness", expand=False))
    stats = calculate_statistics(assignment)
    print_statistics(stats)
    if unassigned_info:
        print_detailed_unassigned_report(unassigned_info)
    else:
        console.print(Panel("[green]All workers assigned successfully![/green]", title="Unassigned Info", expand=False))

print_final_result(final_assignment, final_score, unassigned_info)
