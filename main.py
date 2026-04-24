# main.py
from geometry import GeometryBuilder
from gui import App
from io_utils import parse_rays


def main():

    print("Enter fan rays in format: (1, 0), (0, 1), (-1, -1)")
    #raw = input("Rays: ")
    raw = "(0, 1), (1, 2), (1, 1), (1, 0), (0, -1), (-1, -1)"

    fan = parse_rays(raw)

    print("\nParsed fan:")
    for r in fan:
        print(r)

    geometry = GeometryBuilder(fan).build()

    alpha_dim = geometry.divisor_dimension
    beta_dim = geometry.divisor_dimension

    app = App(resolution=200)

    app.set_geometry(geometry)

    app.run()


if __name__ == "__main__":
    main()