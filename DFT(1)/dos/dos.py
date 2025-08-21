import os
import json
import numpy as np
import matplotlib.pyplot as plt

# 设置基础路径（请根据你的实际路径修改）
base_path = r"C:\Users\Administrator\Desktop\data\DFT\dos"
poscar_path = os.path.join(base_path, "POSCAR")
doscar_path = os.path.join(base_path, "DOSCAR")
outcar_path = os.path.join(base_path, "OUTCAR")  # 可选
output_dir = os.path.join(base_path, "DOS_output")

os.makedirs(output_dir, exist_ok=True)

# 解析 POSCAR 获取元素类型
def parse_poscar(poscar_path):
    with open(poscar_path, 'r') as f:
        lines = f.readlines()
    elements = lines[5].split()
    counts = list(map(int, lines[6].split()))
    atom_types = []
    for elem, count in zip(elements, counts):
        atom_types.extend([elem] * count)
    return atom_types

# 解析 DOSCAR 和 OUTCAR
def parse_doscar(doscar_path, outcar_path):
    if not os.path.exists(doscar_path):
        print(f"Error: {doscar_path} not found.")
        exit(1)

    with open(doscar_path, 'r') as f:
        lines = f.readlines()

    natom = int(lines[0].split()[0])
    efermi = 0.0
    nspin = 1
    lorbit = 0

    if os.path.exists(outcar_path):
        with open(outcar_path, 'r') as f:
            out_lines = f.readlines()
        for line in out_lines:
            if 'E-fermi' in line:
                efermi = float(line.split()[2])
            if 'ISPIN' in line:
                nspin = int(line.split()[2])
            if 'LORBIT' in line:
                lorbit = int(line.split()[2])

    print(f"Number of atoms: {natom}")
    print(f"Fermi energy: {efermi} eV")
    print(f"Spin polarized: {nspin == 2}")
    print(f"LORBIT = {lorbit}")

    start_line = 6
    nl = int(lines[5].split()[2])
    dos_data = []
    for i in range(natom + 1):
        start = start_line + i * (nl + 1)
        end = start + nl
        dos = np.loadtxt(lines[start:end])
        dos_data.append(dos)

    return {
        'natom': natom,
        'efermi': efermi,
        'nspin': nspin,
        'lorbit': lorbit,
        'nl': nl,
        'dos_data': dos_data
    }

# 按元素分类并合并 s/p/d 轨道态密度
def aggregate_element_orbital_pdos(data, atom_types):
    dos_data = data['dos_data']
    nspin = data['nspin']
    efermi = data['efermi']
    lorbit = data['lorbit']

    if lorbit != 11:
        raise ValueError("This script only supports LORBIT=11 for orbital-resolved DOS.")

    element_orbital_pdos = {}

    for i in range(1, data['natom'] + 1):
        atom_dos = dos_data[i]
        energy = atom_dos[:, 0] - efermi

        if nspin == 2:
            s_up = atom_dos[:, 1]
            s_down = atom_dos[:, 2]
            px_up = atom_dos[:, 3]
            px_down = atom_dos[:, 4]
            py_up = atom_dos[:, 5]
            py_down = atom_dos[:, 6]
            pz_up = atom_dos[:, 7]
            pz_down = atom_dos[:, 8]
            dxy_up = atom_dos[:, 9]
            dxy_down = atom_dos[:, 10]
            dyz_up = atom_dos[:, 11]
            dyz_down = atom_dos[:, 12]
            dz2_up = atom_dos[:, 13]
            dz2_down = atom_dos[:, 14]
            dxz_up = atom_dos[:, 15]
            dxz_down = atom_dos[:, 16]
            dx2_y2_up = atom_dos[:, 17]
            dx2_y2_down = atom_dos[:, 18]

            p_up = px_up + py_up + pz_up
            p_down = px_down + py_down + pz_down
            d_up = dxy_up + dyz_up + dz2_up + dxz_up + dx2_y2_up
            d_down = dxy_down + dyz_down + dz2_down + dxz_down + dx2_y2_down
        else:
            s_up = atom_dos[:, 1]
            s_down = np.zeros_like(s_up)
            p_up = atom_dos[:, 2] + atom_dos[:, 3] + atom_dos[:, 4]
            p_down = np.zeros_like(p_up)
            d_up = atom_dos[:, 5] + atom_dos[:, 6] + atom_dos[:, 7] + atom_dos[:, 8] + atom_dos[:, 9]
            d_down = np.zeros_like(d_up)

        elem = atom_types[i - 1]

        if elem not in element_orbital_pdos:
            element_orbital_pdos[elem] = {
                'energy': energy,
                's_up': s_up,
                's_down': s_down,
                'p_up': p_up,
                'p_down': p_down,
                'd_up': d_up,
                'd_down': d_down
            }
        else:
            element_orbital_pdos[elem]['s_up'] += s_up
            element_orbital_pdos[elem]['s_down'] += s_down
            element_orbital_pdos[elem]['p_up'] += p_up
            element_orbital_pdos[elem]['p_down'] += p_down
            element_orbital_pdos[elem]['d_up'] += d_up
            element_orbital_pdos[elem]['d_down'] += d_down

    return element_orbital_pdos

# 保存按元素分类的轨道态密度
def save_element_orbital_pdos(element_orbital_pdos, output_dir):
    for elem, pdos in element_orbital_pdos.items():
        data = {
            "energy": pdos['energy'].tolist(),
            "s_up": pdos['s_up'].tolist(),
            "s_down": (-pdos['s_down']).tolist(),  # 保持与原始代码一致的负值处理
            "p_up": pdos['p_up'].tolist(),
            "p_down": (-pdos['p_down']).tolist(),
            "d_up": pdos['d_up'].tolist(),
            "d_down": (-pdos['d_down']).tolist()
        }

        filename = os.path.join(output_dir, f"{elem}_dos.json")
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)

# 绘制 DOS 图（按轨道和自旋）
def plot_dos_by_orbital(element_orbital_pdos, output_dir):
    fig, axes = plt.subplots(len(element_orbital_pdos), 1, figsize=(8, 4 * len(element_orbital_pdos)))
    if len(element_orbital_pdos) == 1:
        axes = [axes]

    for ax, (elem, pdos) in zip(axes, element_orbital_pdos.items()):
        energy = pdos['energy']
        s_up = pdos['s_up']
        s_down = pdos['s_down']
        p_up = pdos['p_up']
        p_down = pdos['p_down']
        d_up = pdos['d_up']
        d_down = pdos['d_down']

        ax.axhline(0, color='black', lw=1)
        ax.axvline(0, color='red', linestyle='--', alpha=0.5, label='Fermi Level')
        ax.plot(energy, s_up, label='s Spin Up', color='blue')
        ax.plot(energy, -s_down, label='s Spin Down', color='lightblue')
        ax.plot(energy, p_up, label='p Spin Up', color='green')
        ax.plot(energy, -p_down, label='p Spin Down', color='lightgreen')
        ax.plot(energy, d_up, label='d Spin Up', color='orange')
        ax.plot(energy, -d_down, label='d Spin Down', color='gold')
        ax.set_xlabel('Energy - E_F (eV)')
        ax.set_ylabel('DOS (states/eV)')
        ax.set_title(f'{elem} Orbital-resolved DOS')
        ax.legend()
        ax.grid(True)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, 'DOS_orbital_plot.png')
    plt.savefig(plot_path, dpi=300)
    print(f"Plot saved as: {plot_path}")
    plt.show()

# 主程序
def main():
    print("Parsing POSCAR...")
    atom_types = parse_poscar(poscar_path)

    print("Parsing DOSCAR...")
    data = parse_doscar(doscar_path, outcar_path)

    print("Aggregating element orbital PDOS...")
    element_orbital_pdos = aggregate_element_orbital_pdos(data, atom_types)

    print("Saving element orbital PDOS files...")
    save_element_orbital_pdos(element_orbital_pdos, output_dir)

    print("Plotting DOS by orbital...")
    plot_dos_by_orbital(element_orbital_pdos, output_dir)

if __name__ == "__main__":
    main()