def j_overlap(x1, x2, x3, x4):
    """
    Check overlapping of intervals of two elements.
    First element has endpoints x1, x2; second has x3, x4.

    Returns:
        0 if intervals do NOT overlap (no contact)
        1 if intervals overlap (contact)
    """
    # Swap coordinates so that they are always in ascending order
    # for the first element
    if x1 > x2:
        p = x1
        x1 = x2
        x2 = p

    # for the second element
    if x3 > x4:
        p = x3
        x3 = x4
        x4 = p

    # Check interval overlapping

    if x2 == x3 and x1 < x4:
        c = 0  #
    elif x4 == x1 and x3 < x2:
        c = 0  #
    elif x2 < x3 and x3 < x4:
        c = 0  #
    elif x4 < x1 and x1 < x2:
        c = 0  #
    elif x4 < x1 and x1 < x2 and x3 == x4:
        c = 0  #
    elif x2 < x3 and x3 == x4:
        c = 0  #
    elif x4 < x1 and x1 == x2:
        c = 0  #
    elif x1 < x3 and x3 < x4 and x1 == x2:
        c = 0  #
    else:
        c = 1  #
    return c