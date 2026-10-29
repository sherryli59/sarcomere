"""Stable reason labels when resumed files retain historical /meta schemas."""
DEFAULTS = {"0":"not_a_detachment", "1":"geometric_contour", "2":"stochastic",
            "3":"barbed_tip_capsule_penetration", "4":"actn2_myosin_collision",
            "5":"mature_connector_exclusion"}


def reason_names(attrs):
    raw = attrs.get("aa_local_link_detachment_reason_codes", "")
    if isinstance(raw, bytes): raw = raw.decode()
    result = dict(DEFAULTS)
    result.update(item.split("=",1) for item in raw.split(",") if "=" in item)
    return result
