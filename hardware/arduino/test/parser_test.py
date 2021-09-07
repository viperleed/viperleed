import configparser

def main():
    config = configparser.ConfigParser(comment_prefixes='/', allow_no_value=True)
    coeff = '(0, 1)'
    config.read('test_config.ini')
    config.set('energy_calibration', 'coefficients', coeff)
    # config.add_section('Section1')
    # config.set('Section1', 'an_int', '15')
    # config.set('Section1', 'a_bool', 'true')
    # config.set('Section1', 'a_float', '3.1415')
    # config.set('Section1', 'baz', 'fun')
    # config.set('Section1', 'bar', 'Python')
    # config.set('Section1', 'foo', '%(bar)s is %(baz)s!')
    # text = config.get('energy_calibration','coefficients')
    # print(text)
    with open('test_config.ini', 'w') as configfile:
        config.write(configfile)

if __name__ == '__main__':
    main()
