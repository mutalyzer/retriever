from schema import Schema, And, Optional

feature_point = Schema({'type': 'point',
                        'position': And(int, lambda n: n >= 0),
                        Optional('uncertain'): And(bool, True)})

location = Schema({'type': str,
                   'start': feature_point,
                   'end': feature_point,
                   Optional('strand'): int})

feature = Schema({'type': And(str, len),
                  Optional('id'): str,
                  Optional('features'): list,
                  Optional('location'): location,
                  Optional('qualifiers'): dict})


def validate(current):
    feature.validate(current)
    if current.get('features'):
        for sub_feature in current['features']:
            validate(sub_feature)
