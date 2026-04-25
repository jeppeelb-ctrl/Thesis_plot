# main.py
from geometry import GeometryBuilder
from gui import App
from io_utils import parse_rays


def main():

    print("Enter fan rays in format: (1, 0), (0, 1), (-1, -1)")
    #raw = input("Rays: ")
    raw = "(0, 1), (1, 0), (-1, -1), (-1, 0)"

    fan = parse_rays(raw)

    print("\nParsed fan:")
    for r in fan:
        print(r)

    geometry = GeometryBuilder(fan).build()

    app = App()

    app.set_geometry(geometry)

    app.run()


if __name__ == "__main__":
    main()