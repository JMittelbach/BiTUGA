#!/usr/bin/env python3
from typing import Dict, Iterable, List


DEFAULT_PALETTE = ["#c0392b", "#2471a3", "#8e44ad", "#16a085", "#d35400", "#2c3e50"]
MINT_PALETTE = ["#7fb3c7", "#a5cfa1", "#c0d7e3", "#cde7c7", "#aed6c7", "#c6d5cb"]

BASE_COLORS = {
    "female": "#c0392b",
    "male": "#2471a3",
    "trait1": "#c0392b",
    "trait2": "#2471a3",
    "tie": "#7f8c8d",
    "ties": "#7f8c8d",
}


def get_color_map(traits: Iterable[str], palette: str = "default") -> Dict[str, str]:
    """Return a stable color map for the given traits."""
    traits_list: List[str] = list(traits)
    colors = dict(BASE_COLORS)
    pal = MINT_PALETTE if palette == "mint" else DEFAULT_PALETTE
    idx = 0
    for t in traits_list:
        if t in colors:
            continue
        colors[t] = pal[idx % len(pal)]
        idx += 1
    return colors
