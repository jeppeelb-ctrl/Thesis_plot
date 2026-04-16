import random
import time

def blowup_operation(Cone, rays_of_fan, cones_of_fan):
    '''
    This function performs the star subdivision of a cone

    Cone: rays of the cone to be subdivided
    rays_of_fan: the initial rays of the fan
    cones_of_fan: the initial cones of the fan

    Returns: 'new_cone1' and 'new_cone2' are the new cones coming from the subdivison
             'cone_index' is the index of the initial cone, such that we insert
             the new cones in the correct place
    '''
    temp = []
    for i in Cone:
        temp.append(i)
    new_ray = (temp[0][0] + temp[1][0], temp[0][1] + temp[1][1])
    for i in range(len(rays_of_fan)):
        if tuple(temp[0]) == rays_of_fan[i]:
            rays_of_fan.insert(i + 1, tuple(new_ray))

    new_cone1 = [tuple(temp[0]), tuple(new_ray)]
    new_cone2 = [tuple(new_ray), tuple(temp[1])]

    cone_index = cones_of_fan.index(Cone)

    return new_cone1, new_cone2, cone_index


def define_cones(list_of_rays):
    '''
    This function defines the cones of the fan from the rays of the fan

    list_of_rays: is a list containing the vector values of the rays of the fan

    Returns: a list consisting of the cones of the fan
    '''
    cones_of_fan = []

    for i in range(len(list_of_rays)):
        if i == len(list_of_rays) - 1:
            cones_of_fan.append(Cone([list_of_rays[i], list_of_rays[0]]))
        else:
            cones_of_fan.append(Cone([list_of_rays[i], list_of_rays[i + 1]]))

    return cones_of_fan


def blowup_cone(cone_to_blowup, rays_of_fan, cones_of_fan):
    '''
    This function initiates the functions to perform the blow-up of a fan
    '''
    copy_rays = rays_of_fan.copy()
    copy_cones = cones_of_fan.copy()
    new_cone1, new_cone2, cone_index = blowup_operation(cone_to_blowup, copy_rays, copy_cones)
    copy_cones[cone_index] = new_cone1
    copy_cones.insert(cone_index + 1, new_cone2)
    return copy_rays, copy_cones


def hashable_rays(rays):
    return tuple(tuple(ray) for ray in rays)


def is_rays_cached(rays):
    return hashable_rays(rays) in unique_rays


def is_not_useful_rays_cached(rays):
    return hashable_rays(rays) in not_useful_rays


def save_rays_to_cache(rays):
    unique_rays.add(hashable_rays(rays))


def save_not_useful_rays_to_cache(rays):
    not_useful_rays.add(hashable_rays(rays))


def convert_to_jupyter_ray_input(rays):
    str_rays = str(rays)
    str_rays = str_rays.replace("[","")
    str_rays = str_rays.replace("]", "")
    str_rays = str_rays.replace("(", "")
    str_rays = str_rays.replace(")", "")
    str_rays = str_rays.replace(" ", "")

    elements = str_rays.split(",")
    result = elements[0]
    for i in range(1, len(elements)):
        if i % 2 == 0:
            result += ";" + elements[i]
        else:
            result += "," + elements[i]

    return result


def compute_intersection_matrix(number_of_rays, rays):
    '''
    This function computes the intersection matrix of a fan
        and checks if it contains any self-intersection of a specified value

    number_of_rays: the number of rays of the fan
    rays: the rays of the fan

    Returns: 'intersection_matrix' is the intersection matrix of the fan
             'bool_tmp' True if the intersection matrix contains a divisor
             with self-intersection lower than what is desired, and False
             if the intersection matrix does not contain a divisor with self-
             intersection lower than what is desired
    '''
    bool_tmp = False
    intersection_matrix = [[0 for i in range(number_of_rays)] for j in range(number_of_rays)]
    for i in range(number_of_rays):
        if i == number_of_rays - 1:
            left_ray = rays[i - 1]
            right_ray = rays[0]
        else:
            left_ray = list(rays[i - 1])
            right_ray = list(rays[i + 1])
        intersection_number = 0
        for number in range(-15, 15):
            added_rays = [left_ray[0] + right_ray[0], left_ray[1] + right_ray[1]]
            numb_rays = [number * rays[i][0], number * rays[i][1]]
            if added_rays == numb_rays:
                intersection_number = -number
        intersection_matrix[i][i] = intersection_number
        if intersection_number <= -3:
            bool_tmp = True
            break
    return intersection_matrix, bool_tmp

# The initial rays of the fan
ray_input = [(0,1), (1,0), (0,-1), (-1,2)]
cones = define_cones(ray_input)

unique_rays = set()
not_useful_rays = set()

new_rays = []
new_cones = []

# The following randomly performs the star subdivision of the initial fan
# the desired amount of times and then saves the rays of the fan if it
# satisfies the self-intersection criteria, and discards the ones not
# satisfying the criteria

for i in range(500_000):
    current_rays = ray_input
    current_cones = cones

    for _ in range(5):
        r = random.randint(0, len(current_rays) - 1)
        current_rays, current_cones = blowup_cone(current_cones[r], current_rays, current_cones)
    current_bool = False
    if not is_rays_cached(current_rays) and not is_not_useful_rays_cached(current_rays):
        intersection_matrix, current_bool = compute_intersection_matrix(len(current_rays), current_rays)
    if i % 100_000 == 0:
        print(f"iteration {i}: saved rays = {len(new_rays)}, saved -3 rays = {len(not_useful_rays)}")
        print(f"total saved = {len(new_rays) + len(not_useful_rays)}")
    if not current_bool:
        if not is_rays_cached(current_rays) and not is_not_useful_rays_cached(current_rays):
            save_rays_to_cache(current_rays)
            new_rays.append(current_rays)
        else:
            continue
    else:
        if not is_rays_cached(current_rays) and not is_not_useful_rays_cached(current_rays):
            save_not_useful_rays_to_cache(current_rays)

input_rays = []
for k in range(len(new_rays)):
    converted = convert_to_jupyter_ray_input(new_rays[k])
    input_rays.append(converted)