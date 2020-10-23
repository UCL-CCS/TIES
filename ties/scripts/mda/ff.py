import parmed as pmd
import json


def check(top, name):
    # find the atom
    for a in top.atoms:
        if a.name == name:
            atom = a
            break
    print("Atom", a)
    # find all dihedrals that invole our atom
    # print('Angles', len(atom.angles))
    # print('Dihedrals', len(atom.dihedrals))
    return atom


def remove_app_dis_terms(top, app, dis):
    # remove the angles and dihedrals that use atoms coming from both,
    # app and dis sources
    removed_angles = []
    for angle in top.angles[::-1]:
        # if an angle involves an atom fro app and dis, then remove this angle
        atom_names = set([angle.atom1.name, angle.atom2.name, angle.atom3.name])
        if len(atom_names.intersection(app)) > 0 and len(atom_names.intersection(dis)) > 0:
            print('Removing "crossing" angle: ', angle)
            top.angles.remove(angle)
            removed_angles.append(angle)

    removed_dihedrals = []
    for dihedral in top.dihedrals[::-1]:
        # if an angle involves an atom fro app and dis, then remove this angle
        atom_names = set([dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name])
        if len(atom_names.intersection(app)) > 0 and len(atom_names.intersection(dis)) > 0:
            print('Removing "crossing" dihedral: ', dihedral)
            top.dihedrals.remove(dihedral)
            removed_dihedrals.append(dihedral)

    return removed_angles, removed_dihedrals


def uncommon_angles(ta, ma, app, dis):
    print('\n\nAngles not common across the molecules:')
    # take a shallow copy of each
    ta_angles = ta.angles[:]
    ma_angles = ma.angles[:]

    # remove the angles again if they are cross molecules
    for angle in ma_angles[::-1]:
        atom_names = set([angle.atom1.name, angle.atom2.name, angle.atom3.name])
        if len(atom_names.intersection(app)) > 0 and len(atom_names.intersection(dis)) > 0:
            ma_angles.remove(angle)

    not_common = []
    for ta_angle in ta_angles[::-1]:
        for ma_angle in ma_angles[::-1]:
            # if ma_angle in rm_angles:
            # 	break
            one = [ta_angle.atom1.name, ta_angle.atom2.name, ta_angle.atom3.name]
            two = [ma_angle.atom1.name, ma_angle.atom2.name, ma_angle.atom3.name]
            two_r = [ma_angle.atom3.name, ma_angle.atom2.name, ma_angle.atom1.name]

            if one == two or one == two_r:
                ta_angles.remove(ta_angle)
                ma_angles.remove(ma_angle)
                break
    print("left", ta_angles)
    print("right", ma_angles)
    return not_common


def extract_names_di(dihedrals):
    # just for printing
    names = []
    for d in dihedrals:
        names.append([d.atom1.name, d.atom2.name, d.atom3.name, d.atom4.name])
    return names


def uncommon_dihedrals(ta, ma, app, dis):
    print('\n\nCheck which dihedrals are not common')
    # take a shallow copy of each
    ta_dihedrals = ta.dihedrals[:]
    ma_dihedrals = ma.dihedrals[:]

    # remove the angles again if they are cross molecules
    for dihedral in ma_dihedrals[::-1]:
        atom_names = set([dihedral.atom1.name, dihedral.atom2.name, dihedral.atom3.name, dihedral.atom4.name])
        if len(atom_names.intersection(app)) > 0 and len(atom_names.intersection(dis)) > 0:
            ma_dihedrals.remove(dihedral)

    not_common = []
    for ta_angle in ta_dihedrals[::-1]:
        for ma_angle in ma_dihedrals[::-1]:
            one = [ta_angle.atom1.name, ta_angle.atom2.name, ta_angle.atom3.name, ta_angle.atom4.name]
            two = [ma_angle.atom1.name, ma_angle.atom2.name, ma_angle.atom3.name, ma_angle.atom4.name]
            two_r = [ma_angle.atom4.name, ma_angle.atom3.name, ma_angle.atom2.name, ma_angle.atom1.name]

            if one == two or one == two_r:
                ta_dihedrals.remove(ta_angle)
                ma_dihedrals.remove(ma_angle)
                break
    print("Left:\n", '\n'.join('-'.join(d) for d in extract_names_di(ta_dihedrals)))
    print("\nRight:\n", '\n'.join('-'.join(d) for d in extract_names_di(ma_dihedrals)))
    return not_common


# load the json
with open('joint_meta_fep.json') as f:
    suptop = json.load(f)

    app = set(suptop['appearing'])
    print('app', app)
    dis = set(suptop['disappearing'])
    print('dis', dis)

# left C7 is connected to CL8
l_aname = 'C7'
lt = pmd.load_file('left_top/sys.top')

# C59 is connected to BR1
r_aname = 'C59'
rt = pmd.load_file('right_top/sys.top')

mt = pmd.load_file('morph_top/sys.top')
print('mt angles', len(mt.angles))
rm_angles, rm_dihedrals = remove_app_dis_terms(mt, app, dis)
print('mt after angles', len(mt.angles))
print('\n\n')


# compare mt to lt and check if there is anything extra besides LT
# these extra stuff in mt-lt comparison should be touching on the rt only
def in_one_oranother(mt, alonet, app):
    for mtang in mt.angles:
        mt_ang_names = [mtang.atom1.name, mtang.atom2.name, mtang.atom3.name]
        found_in_alonet = False
        for alonetang in alonet.angles:
            alonet_ang_names = [alonetang.atom1.name, alonetang.atom2.name, alonetang.atom3.name]
            if mt_ang_names == alonet_ang_names or mt_ang_names == alonet_ang_names[::-1]:
                found_in_alonet = True

        if not found_in_alonet:
            # app is rt, check if it there
            found_in_rt = False
            for aname in mt_ang_names:
                if aname in app:
                    found_in_rt = True

            if not found_in_rt:
                print('mt ang not in one or dis/app', mt_ang_names)


print('mt lt app')
in_one_oranother(mt, lt, app)
print('mt rt dis')
in_one_oranother(mt, rt, dis)

# than do the same for rt, any extra stuff

# check the standalone
# lt_c7 = check(lt, l_aname)
# rt_c59 = check(rt, r_aname)
# # check the morph
# mt_c7 = check(mt, l_aname)
# mt_c59 = check(mt, r_aname)
# # which angles are not in common?
# uncommon_angles(lt_c7, mt_c7, app, dis)
# uncommon_dihedrals(lt_c7, mt_c7, app, dis)
# print('\n\nangles')
# print('lt_c7', lt_c7.angles)
# print('mt_c7', mt_c7.angles)
# print('mt angles at the end', len(mt.angles))
