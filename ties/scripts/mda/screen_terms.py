import os
import json

import parmed as pmd

# load the json
with open('joint_meta_fep.json') as f:
    suptop = json.load(f)

    app = set(suptop['appearing'])
    print('app', app)
    dis = set(suptop['disappearing'])
    print('dis', dis)

# load the morph file
ligtop = pmd.load_file('lig/sys_solv.top')
print ('all angles', len(ligtop.angles))

for angle in ligtop.angles[::-1]:
    # if an angle involves an atom fro app and dis, then remove this angle
    atom_names = set([angle.atom1.name, angle.atom2.name, angle.atom3.name])
    if len(atom_names.intersection(app)) > 0 and len(atom_names.intersection(dis)) > 0:
        print('Removing angle: ', angle)
        ligtop.angles.remove(angle)

for dihedral in ligtop.dihedrals[::-1]:
    # if an angle involves an atom fro app and dis, then remove this angle
    atom_names = set([dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name])
    if len(atom_names.intersection(app)) > 0 and len(atom_names.intersection(dis)) > 0:
        print('Removing dihedral: ', dihedral)
        ligtop.dihedrals.remove(dihedral)

# backup the original file
os.rename('lig/sys_solv.top', 'lig/sys_solv.top.original')
ligtop.save('lig/sys_solv.parm7')
os.rename('lig/sys_solv.parm7', 'lig/sys_solv.top')