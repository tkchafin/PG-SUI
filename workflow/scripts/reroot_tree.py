import argparse
import toytree

def reroot_tree(tree, rooted="out.rooted.tre", outgroup_wildcard="out"):
    t=toytree.tree(tree)
    try:
        rt=t.root(wildcard=outgroup_wildcard)
        rt.write(rooted, tree_format=5)
        return(rt)
    except Exception:
        t.write(rooted, tree_format=5)
        return(None)

def main():
    parser = argparse.ArgumentParser(description='Reroot a tree using an outgroup.')
    parser.add_argument('-i', '--input', type=str, help='Input tree file.', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output tree file.', required=True)
    parser.add_argument('-O', '--outgroup', type=str, default="out", help='Wildcard for the outgroup.')
    
    args = parser.parse_args()
    reroot_tree(args.input, args.output, args.outgroup)

if __name__ == "__main__":
    main()
