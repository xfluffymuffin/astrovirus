import svgutils.compose as sc
import os.path

svg_files = []

for dirpath, dirnames, filenames in os.walk(r"path_to_directory"): # Директория с объединяемыми svg-изображениями
    for filename in [f for f in filenames if f.endswith(".svg")]:
        svg_files.append(os.path.join(dirpath, filename))
        # print(os.path.join(dirpath, filename))


# Размер холста (в дюймах, по умолчанию SVG единицы — pt)
canvas_width = 12
canvas_height = 19

layout = sc.Figure(
    canvas_width * 72,  # Перевод в pt (1 дюйм = 72 точки)
    canvas_height * 72,
    sc.Panel(sc.SVG(svg_files[0]).move(0, 0)),
    sc.Panel(sc.SVG(svg_files[1]).move(420, 0)),

    sc.Panel(sc.SVG(svg_files[2]).move(0, 340)),
    sc.Panel(sc.SVG(svg_files[3]).move(420, 340)),

    sc.Panel(sc.SVG(svg_files[4]).move(0, 2 * 340)),
    sc.Panel(sc.SVG(svg_files[5]).move(420, 2 * 340)),

    sc.Panel(sc.SVG(svg_files[6]).move(210, 3 * 340)))

    #sc.Panel(sc.SVG(svg_files[6]).move(0, 3 * 340)),
    #sc.Panel(sc.SVG(svg_files[7]).move(420, 3 * 340)))

layout.save("filename.svg")
