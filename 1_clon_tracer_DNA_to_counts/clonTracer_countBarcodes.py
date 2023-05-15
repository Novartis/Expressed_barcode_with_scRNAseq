__author__ = 'krishvi7'
import sys
import getopt
from modules import classes

clt_object = classes.ClonTracerCountBarcodes()


def main(argv):
    arg_single = ":".join(clt_object.args_dict.values())+":" + "".join(clt_object.args_dict_paramless.values())
    arg_verbose = map(lambda x: x + "=", clt_object.args_dict.keys())
    arg_verbose.extend(map(lambda x: x, clt_object.args_dict_paramless.keys()))

    try:
        opts, args = getopt.getopt(argv, arg_single, arg_verbose)
    except getopt.GetoptError as exp:
        classes.Common.print_error("Invalid input argument: "+exp.msg)

    if len(opts) == 0:
        clt_object.usage()

    for opt, arg in opts:
        if opt in map(lambda x: "-"+x, clt_object.args_dict.values()) or \
           opt in map(lambda x: "--"+x, clt_object.args_dict.keys()):
            if opt in map(lambda x: "-"+x, clt_object.args_dict.values()):
                clt_member = clt_object.args_dict.keys()[clt_object.args_dict.values().index(opt[1:])]
            else:
                clt_member = opt[2:]

            try:
               tryout = classes.Common().dataType_dict[clt_member](arg)
            except ValueError as val_err:
                classes.Common.print_error("Invalid argument "+arg+" for parameter "+opt)

            if clt_member in classes.Common().range_dict.keys():
                if classes.Common().dataType_dict[clt_member](arg) not in classes.Common().range_dict[clt_member]:
                    raise classes.Common.print_error("Argument for "+opt+" must be " +
                                                         ",".join(map(str, classes.Common().range_dict[clt_member])))

            setattr(clt_object, clt_member, classes.Common().dataType_dict[clt_member](arg))

        elif opt in map(lambda x: "-"+x, clt_object.args_dict_paramless.values()) or \
           opt in map(lambda x: "--"+x, clt_object.args_dict_paramless.keys()):
            if opt in map(lambda x: "-"+x, clt_object.args_dict_paramless.values()):
                clt_member = clt_object.args_dict_paramless.keys()[clt_object.args_dict_paramless.values().index(opt[1:])]
            else:
                clt_member = opt[2:]

            setattr(clt_object, clt_member, classes.Common().paramless_dict[clt_member])
        else:
            classes.Common.print_error("Invalid input argument "+opt)

    try:
        clt_object.check_input_files()
        clt_object.count_barcodes()
        clt_object.simplify_barcodes()
    except Exception:
        raise

if __name__ == "__main__":
    main(sys.argv[1:])
