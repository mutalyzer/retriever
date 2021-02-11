def make_location(start, end=None, strand=None):
    if end is not None:
        location = {
            "type": "range",
            "start": {"type": "point", "position": int(start)},
            "end": {"type": "point", "position": int(end)},
        }
    else:
        location = {"type": "point", "position": int(start)}
    if strand is not None:
        location["strand"] = int(strand)
    return location
