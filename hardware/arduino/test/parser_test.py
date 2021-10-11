# for testing only, delete if done testing
import configparser
import ast

def main():
    config = configparser.ConfigParser(comment_prefixes='/', allow_no_value=True)
    coeff = '(0, 1)'
    file_name = 'test_config.ini'
    try:
        f = open(file_name, 'r')
        f.close()
        config.read(file_name)
    except OSError:
        raise FileNotFoundError(f"Could not open/read file: {file_name}. "
                                "Check if this file exists and is in the "
                                "correct folder.")
        return
    # config.set('energy_calibration', 'coefficients', coeff)
    # config.add_section('Section1')
    # config.set('Section1', 'an_int', '15')
    # config.set('Section1', 'a_bool', 'true')
    # config.set('Section1', 'a_float', '3.1415')
    # config.set('Section1', 'baz', 'fun')
    # config.set('Section1', 'bar', 'Python')
    # config.set('Section1', 'foo', '%(bar)s is %(baz)s!')
    # text = config.get('energy_calibration','coefficients')
    # print(text)
    prim_set = ast.literal_eval(config.get('devices','secondary_controllers'))
    for controller in prim_set:
        print(controller[0])
    cam_set = ast.literal_eval(config.get('devices','cameras'))
    for camera in cam_set:
        print(camera)

if __name__ == '__main__':
    main()
