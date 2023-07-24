import os
import glob
from rmsdfn import get_coor_from_sdf, accurate_rmsd

from argparse import FileType, ArgumentParser
parser = ArgumentParser()
parser.add_argument('--input_path', type=str, default='input', help='Path to folder with testing ligand-protein structures')
parser.add_argument('--diffdock_result_path', type=str, default='diffdock/results0129', help='Path to folder with testing ligand-protein structures')

parser.add_argument('--results_path', type=str, default='results', help='Path to save results')

parser.add_argument('--seed', type=str, default='666', help='Random seed for docking')
parser.add_argument('--exhaustiveness', type=str, default='32', help='exhaustiveness for docking searching')
parser.add_argument('--cnn_scoring', type=str, default='rescore', choices=['rescore','none'], help='whether or not to use cnn for rescoring')
parser.add_argument('--whole_ptn', action='store_true', help='perform whole protein docking')
parser.add_argument('--gt_pocket', action='store_true', help='use the ground truth pocket for docking')
parser.add_argument('--ori_pocket', action='store_true', default=False, help='use exact pocket of Diffdock output without normalize')
parser.add_argument('--pocket_size', type=float, default=30, help='normalized pocket size in Angstrom')
parser.add_argument('--lig_num', type=int, help='which ligand index to dock, useful for parallelizing this process across a cluster')
parser.add_argument('--cpus', type=int, help='how many cpus for gnina to use, if not specified then GNINA will detect')

args = parser.parse_args()

input_path = args.input_path

seed = args.seed
if args.cnn_scoring == 'rescore':
    sf = 'gnina'
else:
    sf = 'vina'
# mode = "detail"
mode = args.exhaustiveness
normal_pocket = not args.ori_pocket and not args.gt_pocket and not args.whole_ptn
pocket_size = args.pocket_size
if not os.path.isdir(args.results_path):
    os.mkdir(args.results_path)
if normal_pocket:
    output_path = os.path.join(args.results_path, mode+'_'+sf+f"_diffdock_box_pdbqt_{sf}_seed"+str(seed) +"_normal"+str(pocket_size))
elif args.gt_pocket:
    output_path = os.path.join(args.results_path, mode+'_'+sf+f"_diffdock_box_pdbqt_{sf}_seed"+str(seed) +"_gt_pock")
else:
    output_path = os.path.join(args.results_path, mode+'_'+sf+f"_diffdock_box_pdbqt_{sf}_seed"+str(seed))

if args.whole_ptn:
    output_path = os.path.join(args.results_path, mode+'_'+sf+f"_diffdock_box_pdbqt_{sf}_seed"+str(seed)+ "_whole_protein")

os.system("mkdir "+output_path)

out_csv = os.path.join(args.results_path, "rmsd_results_{}_{}_diffdock_box_pdbqt_{}_seed{}.csv".format(mode,sf,sf,seed))
if normal_pocket:
    out_csv = out_csv.replace('.csv','_normal.csv')
elif args.gt_pocket:
    out_csv = out_csv.replace('.csv','_gt_pock.csv')
elif args.whole_ptn:
    out_csv = out_csv.replace('.csv','_whole_protein.csv')
with open(out_csv,"w") as f:
    f.write("pdb,top1_rmsd,top5_rmsd\n")
    
configs = sorted(glob.glob(os.path.join(args.diffdock_result_path,"*/rank1.sdf")),reverse=True)


if args.lig_num is not None:
    configs = [configs[args.lig_num]]
    print(len(configs))
for config in configs:
    out_ext = '.pdbqt'
    pdb = config.split('/')[-2].split('-')[-2]
    with open(config, "r") as f:
        lines = f.readlines()
    
    ligand = os.path.join(input_path,'{}/{}_ligand.pdbqt'.format(pdb,pdb))    
    ligand_sdf = ligand.replace('.pdbqt','.sdf')
    if os.path.isfile(ligand_sdf):
        ligand = ligand_sdf
        out_ext = '.sdf'
    
    ori_coords, _ = get_coor_from_sdf(config)
    
    center_max = [-1e4,-1e4,-1e4]
    center_min = [1e4,1e4,1e4]
    center = [0,0,0]
    size = [0,0,0]
    for co in ori_coords:
        for i in range(3):
            center_max[i] = max(center_max[i], co[i])
            center_min[i] = min(center_min[i], co[i])
    
    for i in range(3):
        center[i] = (center_max[i] + center_min[i]) / 2
        size[i] = (center_max[i] - center_min[i])  + 5
    
    
    new_config = os.path.join(output_path,pdb+"_config.txt")
    with open(new_config,"w") as f:
        if args.gt_pocket:
            f.write("autobox_ligand = {}\n".format(ligand))
        elif args.whole_ptn:
            f.write("autobox_ligand = {}/{}_protein.pdbqt\n".format(os.path.join(input_path,pdb),pdb))
        else:
            f.write("center_x = {}\ncenter_y = {}\ncenter_z = {}\n".format(center[0],center[1],center[2]))
            if normal_pocket:
                f.write("size_x = {}\nsize_y = {}\nsize_z = {}\n".format(pocket_size,pocket_size,pocket_size))
            else:
                f.write("size_x = {}\nsize_y = {}\nsize_z = {}\n".format(min(size[0],60),min(size[1],60),min(size[2],60)))
        f.write("cnn_scoring = {}\n".format(args.cnn_scoring))
        f.write("exhaustiveness = {}\n".format(mode.split('_')[0]))
        
        f.write("receptor = {}/{}_protein.pdbqt\nligand = {}\n".format(os.path.join(input_path,pdb),pdb,ligand))
        f.write("out = {}\nseed = {}\nnum_modes = 9\nverbosity = 1\n".format(os.path.join(output_path,pdb+f"_ligand_out{out_ext}"), seed))
        if args.cpus is not None:
            f.write("cpu = {}\n".format(args.cpus))
        
    if not os.path.isfile(os.path.join(output_path,pdb+f"_ligand_out{out_ext}")):
        os.system("gnina --config "+new_config)
    if not os.path.isfile(os.path.join(output_path,pdb+f"_ligand_out{out_ext}")):
        print("error in "+ pdb)
        with open(out_csv, "a") as out_f:
            out_f.write(pdb+",999\n")
        continue
    rmsds = accurate_rmsd(os.path.join(output_path, pdb+f"_ligand_out{out_ext}"), ligand)[:5]
    top1_rmsd = rmsds[0]
    top5_rmsd = min(rmsds)
    with open(out_csv, "a") as out_f:
        out_f.write(pdb+','+str(top1_rmsd)+','+str(top5_rmsd)+"\n")
        
