import os
import sys
import argparse
import bte
from pango_aliasor.aliasor import Aliasor
import dendropy

global_aliasor = Aliasor()

def process_mstr(mstr):
    """Read a mutation string and return the chromosome, location, reference, and alternate alleles."""
    if ":" in mstr:
        chro = mstr.split(":")[0]
        data = mstr.split(":")[1]
    else:
        chro = None
        data = mstr
    if data[0].isdigit():
        loc = int(data[:-1])
        ref = None
        alt = data[-1]
    else:
        loc = int(data[1:-1])
        ref = data[0]
        alt = data[-1]
    return chro, loc, ref, alt

def compute_mutation_weight(node, mutweights):
    if len(mutweights) == 0:
        if node.branch_length < 0:
            print(f"WARNING: Negative branch length detected on node {node.id}! Treating as 0...")
            return 0
        return node.branch_length
    dist = 0
    for m in node.mutations:
        _, loc, _, alt = process_mstr(m)
        mweight = max([mutweights.get((loc, alt, None), 0), mutweights.get((loc, alt, node.id), 0)])
        dist += mweight
    return dist

def dists_to_root(node, mutweights={}):
    nodes = {node.id: 0}
    def recursive_dists_to_roots(snode):
        bweight = nodes[snode.id]
        for child in snode.children:
            dist = bweight + compute_mutation_weight(child, mutweights)
            nodes[child.id] = dist
            recursive_dists_to_roots(child)
    recursive_dists_to_roots(node)
    return nodes

def get_sum_and_count(rbfs, ignore=set(), mutweights={}, sampleweights={}):
    sum_and_count_dict = {}
    leaf_count = 0
    for node in rbfs:
        if node.is_leaf():
            leaf_count += 1
            if node.id not in ignore:
                if len(sampleweights) == 0:
                    count = 1
                else:
                    count = float(sampleweights.get(node.id, 0))
                sum_and_count_dict[node.id] = (compute_mutation_weight(node, mutweights), count)
        else:
            total_count = 0
            total_sum = 0
            for child in node.children:
                sumtc = sum_and_count_dict.get(child.id, None)
                if sumtc is None:
                    continue
                total_count += sumtc[1]
                total_sum += sumtc[0]
            if total_count > 0:
                sum_and_count_dict[node.id] = (total_sum + compute_mutation_weight(node, mutweights) * total_count, total_count)
    return sum_and_count_dict, leaf_count

def evaluate_candidate(a, nid, sum_and_counts, dist_to_root, minimum_size=0, minimum_distinction=0):
    node_sum, node_count = sum_and_counts.get(nid, [0, 0])
    if node_count <= minimum_size:
        return 0
    if node_sum == 0 or node_count <= 0:
        return 0
    candidate_to_parent = dist_to_root[nid] - dist_to_root[a]
    if candidate_to_parent < minimum_distinction:
        return 0
    mean_distances = node_sum / node_count
    if (mean_distances + candidate_to_parent) == 0:
        candidate_value = 0
    else:
        candidate_value = node_count * candidate_to_parent / (mean_distances + candidate_to_parent)
    return candidate_value

def evaluate_lineage(t, dist_to_root, anid, candidates, sum_and_count, minimum_size=0, minimum_distinction=0, banned=set()):
    good_candidates = []
    for c in candidates:
        if not c.is_leaf() and c.id not in banned:
            cscore = evaluate_candidate(anid, c.id, sum_and_count, dist_to_root, minimum_size, minimum_distinction)
            if cscore > 0:
                good_candidates.append((cscore, c))
    if len(good_candidates) == 0:
        return (0, None)
    return max(good_candidates, key=lambda x: x[0])

def get_skipset(t, annotes):
    skip = set()
    for lin, nid in annotes.items():
        ancestors = t.rsearch(nid, True)
        for anc in ancestors:
            skip.add(anc.id)
    return skip

def get_outer_annotes(t, annotes):
    outer_annotes = {}
    for l in t.get_leaves():
        mann = l.most_recent_annotation()
        for a in mann:
            if a is not None and a not in outer_annotes and a in annotes:
                outer_annotes[a] = annotes[a]
        if len(outer_annotes) == len(annotes):
            break
    return outer_annotes

def parse_mutweights(mutweights_file):
    mutweights = {}
    with open(mutweights_file) as f:
        for line in f:
            line = line.strip()
            if line == "":
                continue
            if line[0] == "#":
                continue
            parts = line.split()
            _, loc, _, alt = process_mstr(parts[0])
            if len(parts) == 3:
                mutweights[(loc, alt, parts[2])] = float(parts[1])
            else:
                mutweights[(loc, alt, None)] = float(parts[1])
    if len(mutweights) == 0:
        print("ERROR: Mutation weight file indicated found empty!")
        exit(1)
    return mutweights

def build_annotation_network(t, rawann):
    annd = {}
    for ann, nid in rawann.items():
        pnode = t.get_node(nid).parent
        if pnode is not None:
            parents = pnode.most_recent_annotation()
        else:
            parents = []
        for p in parents:
            if ann not in annd:
                annd[ann] = [p]
            else:
                annd[ann].append(p)
    return annd

def read_samples_weights(sfile):
    samples = {}
    with open(sfile) as inf:
        for entry in inf:
            spent = entry.strip().split()
            if len(spent) == 1:
                samples[spent[0]] = 1
            else:
                samples[spent[0]] = spent[1]
    return samples

def filter_annotes(t, annotes, selection):
    filtered = {}
    for ann, nid in annotes.items():
        ancestry = t.rsearch(nid, True)
        for a in ancestry:
            if selection in a.annotations:
                filtered[ann] = nid
    return filtered

def parse_aaweights(aaf):
    aad = {}
    with open(aaf) as inf:
        for entry in inf:
            gene, sitel, weight = entry.strip().split()
            site = int(sitel[:-1])
            state = sitel[-1]
            aad[(gene, site, state)] = float(weight)
    return aad

def argparser():
    parser = argparse.ArgumentParser(description="Propose sublineages for existing lineages based on relative representation concept.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to annotate.')
    parser.add_argument("-c", "--clear", action='store_true', help='Clear all current annotations and apply a level of serial annotations to start with.')
    parser.add_argument("-r", "--recursive", action='store_true', help='Recursively add additional sublineages to proposed lineages.')
    parser.add_argument("-o", "--output", help='Path to output protobuf, if desired.', default=None)
    parser.add_argument("-d", "--dump", help="Print proposed sublineages to a table.", default=None)
    parser.add_argument("-l", "--labels", help="Print lineage and sample associations to a table formatted for matUtils annotate -c.", default=None)
    parser.add_argument("-t", "--distinction", help="Require that lineage proposals have at least t mutations distinguishing them from the parent lineage or root.", type=int, default=1)
    parser.add_argument("-m", "--minsamples", help="Require that each lineage proposal represent at least m total sample weight (without special weighting, the number of samples).", type=int, default=10)
    parser.add_argument("-w", "--mutweights", help="Path to an optional two (or three) column space-delimited containing mutations and weights (and nodes) to use to weight lineage choices.", default=None)
    parser.add_argument("-y", "--aaweights", help="Path to an optional three column space-delimited containing amino acid changes and weights to use to weight lineage choices. Requires --gtf and --reference to be set. Changes not included will be weighted as 1.", default=None)
    parser.add_argument("-g", "--gene", help='Consider only mutations in the indicated gene. Requires that --gtf and --reference be set.', default=None)
    parser.add_argument("-s", "--missense", action='store_true', help="Consider only missense mutations. Requires that --gtf and --reference be set.")
    parser.add_argument("-u", "--cutoff", help="Stop adding serial lineages when at least this proportion of samples are covered. Default 0.95", type=float, default=0.95)
    parser.add_argument("-f", "--floor", help="Minimum score value to report a lineage. Default 0", type=float, default=0)
    parser.add_argument("--gtf", help="Path to a gtf file to apply translation. Use with --reference.")
    parser.add_argument("--reference", help='Path to a reference fasta file to apply translation. Use with --gtf.')
    parser.add_argument("-v", "--verbose", help='Print status updates.', action='store_true')
    parser.add_argument("-a", "--annotation", help='Choose a specific lineage, and its sublineages, to propose new sublineages for.', default=None)
    parser.add_argument("-p", "--samples", help='Path to a space-delimited file containing samples and weights in the first and second columns. If used, samples not included in this file will be ignored.', default=None)
    # parser.add_argument("--label_file", help='Path to the label file to replace tree tip labels.', default=None)  # Add this line
    return parser

def propose(args):
    t = bte.MATree(args.input)
    mutweights = {}
    if args.gene == 'ORF1a' or args.gene == 'ORF1b':
        print("WARNING: ORF1a and ORF1b are treated as a unified ORF1ab for purposes of haplotype identification due to complexities with redundant counting and translation implementation.")
        args.gene = "ORF1ab"
    if args.gtf is not None and args.reference is not None:
        if args.verbose:
            print("Performing tree translation and setting weights for mutations based on amino acid changes.")
        aaweights = {}
        if args.aaweights is not None:
            print("Retrieving amino acid change weightings.")
            aaweights = parse_aaweights(args.aaweights)
        translation = t.translate(fasta_file=args.reference, gtf_file=args.gtf)
        for nid, aav in translation.items():
            for aa in aav:
                if args.missense and aa.is_synonymous():
                    continue
                if args.gene is None or args.gene == "None" or aa.gene == args.gene:
                    mutweights[(int(aa.nt_index), aa.alternative_nt, nid)] = aaweights.get((aa.gene, aa.aa_index, aa.aa), 1)
        if len(mutweights) == 0:
            raise ValueError("No mutations have weights after translation! Check parameters")
    if args.mutweights is not None:
        mutweights.update(parse_mutweights(args.mutweights))
    if args.verbose:
        print("Considering {} mutations to have weight.".format(len(mutweights)))
    if args.dump is not None:
        dumpf = open(args.dump, 'w+')
    if args.clear:
        t.apply_node_annotations({node.id: [] for node in t.depth_first_expansion()})
    try:
        cannotes = t.dump_annotations()
    except:
        cannotes = t.get_annotations()  # replacement function in newer versions of bte
    annotes = {}
    for k, v in cannotes.items():
        annotes[global_aliasor.uncompress(k)] = v
    if args.annotation is not None:
        if args.clear:
            print("ERROR: Cannot select lineages (-a) while clearing lineages (-c)!")
            exit(1)
        if args.verbose:
            print("Finding annotations that are descendants of {}.".format(args.annotation))
        annotes = filter_annotes(t, annotes, args.annotation)
        if args.verbose:
            print("Found {} annotations to check for sublineages.".format(len(annotes)))
    if args.clear:
        assert len(annotes) == 0
    ann_net = build_annotation_network(t, annotes)
    original_annotations = set(annotes.keys())
    global_used_nodes = get_skipset(t, annotes)
    if len(annotes) == 0:
        if args.verbose and not args.clear:
            print("No lineages found in tree; starting from root.")
        annotes = {'L': t.root.id}
    else:
        if args.verbose:
            print("{} annotations found in the tree; identifying candidates for subdivision.".format(len(annotes)))
        annotes = get_outer_annotes(t, annotes)
        if args.verbose:
            print("{} outer annotations found in the tree; identifying sublineages.".format(len(annotes)))
    if args.verbose:
        print("Tree contains {} annotated lineages initially ({} nodes disregarded to prevent retroactive parent assignment).".format(len(annotes), len(global_used_nodes)))
    if args.dump is not None:
        print("parent\tparent_nid\tproposed_sublineage\tproposed_sublineage_nid\tproposed_sublineage_score\tproposed_sublineage_size", file=dumpf)
    outer_annotes = annotes
    global_labeled = set()
    sample_weights = {}
    if args.samples is not None:
        allsamples = t.get_leaves_ids()
        sample_weights = read_samples_weights(args.samples)
        for s in allsamples:
            if s not in sample_weights:
                global_labeled.add(s)
        if args.verbose:
            print("{} samples given weights; ignoring {} samples".format(len(sample_weights), len(global_labeled)))
    level = 1
    while True:
        if args.verbose:
            print("Level: ", level)
        new_annotes = {}
        used_nodes = global_used_nodes.copy()
        for ann, nid in outer_annotes.items():
            serial = 1
            rbfs = t.breadth_first_expansion(nid, True)  # takes the name
            if len(sample_weights) == 0:
                parent_leaf_count = len([n for n in rbfs if n.is_leaf()])
            else:
                parent_leaf_count = len([n for n in rbfs if n.id in sample_weights])
            if parent_leaf_count == 0:
                if args.verbose:
                    print("No samples descended from {} have weight; continuing".format(ann))
                continue
            current_child_lineages = {k: v for k, v in annotes.items() if ann in ann_net.get(k, [])}
            labeled = global_labeled.copy()
            for lin, cnid in current_child_lineages.items():
                for s in t.get_leaves_ids(cnid):
                    labeled.add(s)
            if len(current_child_lineages) > 0 and args.verbose:
                print("Found {} child lineages preexisting for lineage {}; {} samples prelabeled from {} total ({}%)".format(len(current_child_lineages), ann, len(labeled) - len(global_labeled), parent_leaf_count, 100 * (len(labeled) - len(global_labeled)) / parent_leaf_count))
            dist_root = dists_to_root(t.get_node(nid), mutweights)  # needs the node object, not just the name
            while True:
                scdict, leaf_count = get_sum_and_count(rbfs, ignore=labeled, mutweights=mutweights, sampleweights=sample_weights)
                best_score, best_node = evaluate_lineage(t, dist_root, nid, rbfs, scdict, args.minsamples, args.distinction, used_nodes)
                if best_score <= args.floor:
                    break
                if ann[:5] == 'auto.':
                    prefix = ann
                else:
                    prefix = "auto." + ann
                newname = prefix + "." + str(serial)
                while newname in original_annotations or newname.lstrip("auto.") in original_annotations:
                    serial += 1
                    newname = prefix + '.' + str(serial)
                for anc in t.rsearch(best_node.id, True):
                    used_nodes.add(anc.id)
                new_annotes[newname] = best_node.id
                leaves = t.get_leaves_ids(best_node.id)
                if args.dump is not None:
                    print("{}\t{}\t{}\t{}\t{}\t{}".format(ann, nid, newname, best_node.id, str(best_score), len(leaves)), file=dumpf)
                for l in leaves:
                    labeled.add(l)
                if len(labeled) >= leaf_count * args.cutoff:
                    break
                serial += 1
                if args.verbose:
                    print("Annotated lineage {} as descendent of {} from level {} with {} descendents".format(newname, ann, level, len(leaves)))
        if not args.recursive:
            annotes.update(new_annotes)
            break
        elif len(new_annotes) == 0:
            break
        else:
            annotes.update(new_annotes)
            outer_annotes = new_annotes
            level += 1
    if args.verbose:
        print("After sublineage annotation, tree contains {} annotated lineages.".format(len(annotes)), file=sys.stderr)

    if args.output is not None:
        annd = {}
        for k, v in annotes.items():
            try:
                k = global_aliasor.compress(k)
            except:
                pass
            if v not in annd:
                annd[v] = []
            if len(annd[v]) == 2:
                annd[v][1] = k
            else:
                annd[v].append(k)
        t.apply_node_annotations(annd)
        labels = {}
        for ann, nid in annotes.items():
            try:
                ann = global_aliasor.compress(ann)
            except:
                pass
            for leaf in t.get_leaves_ids(nid):
                if leaf not in labels:
                    labels[leaf] = [ann]
                else:
                    labels[leaf].append(ann)
            # with open(args.labels, 'w+') as f:
            #     for l, v in labels.items():
            #         for ann in v:
            #             print("{}\t{}".format(ann, l), file=f)
        replace_labels(t, labels)  # Add this line
        t.save_pb(args.output)
    if args.dump is not None:
        dumpf.close()
    if args.labels is not None:
    #     labels = {}
    #     for ann, nid in annotes.items():
    #         try:
    #             ann = global_aliasor.compress(ann)
    #         except:
    #             pass
    #         for leaf in t.get_leaves_ids(nid):
    #             if leaf not in labels:
    #                 labels[leaf] = [ann]
    #             else:
    #                 labels[leaf].append(ann)
        with open(args.labels, 'w+') as f:
            for l, v in labels.items():
                for ann in v:
                    print("{}\t{}".format(ann, l), file=f)

def replace_labels(tree, labels):
    # with open(label_file, 'r') as f:
    #     labels = {}
    #     for line in f:
    #         lineage, tip = line.strip().split()
    #         if lineage not in labels:
    #             labels[lineage] = []
    #         labels[lineage].append(tip)
    
    print(labels)
    label_counts = {} 
    
    for key in labels:
        base_label = labels[key][0]  
        if base_label not in label_counts:
            label_counts[base_label] = 1  
        else:
            label_counts[base_label] += 1  
        
        labels[key] = [f"{base_label}_{label_counts[base_label]}"]
    print(labels)
    tree.apply_node_annotations(labels)
    # print(tree.get_leaves_ids())

def main():
    parser = argparser()
    args = parser.parse_args()
    propose(args)

if __name__ == "__main__":
    main()
