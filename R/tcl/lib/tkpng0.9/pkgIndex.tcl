#
# Tcl package index file
#

proc loadTkPNGDll {dir plat} {
    set opwd [pwd]
    if {[regexp Windows $plat]} {
        cd  $dir
        # package require Tk
        #tk_messageBox -type ok -message "Trying dll load ...!"
        set success false
        catch {
            load tkpng-Windows-amd64.dll tkpng
            #tk_messageBox -type ok -message "Yeah tkpng-Windows-amd64.dll loaded!!"
            set success true
        } 
        if {!$success} {
            catch {
                load tkpng-Windows-x86.dll tkpng
                #tk_messageBox -type ok -message "Yeah tkpng-Windows-x86.dll loaded!!"
                set success true
            }
        }
        #tk_messageBox -type ok -message "Success?? Success = $success !!"
        cd $opwd
    } else {
        load [file join $dir tkpng-${plat}[info sharedlibextension]] tkpng
    }
    package provide tkpng 0.9
}
if {[catch {package require Tcl 8.2}]} return
set plat [lindex $::tcl_platform(os) 0]-$::tcl_platform(machine)
if {$::tcl_platform(os) eq "Darwin"} {
    # untested
    set plat Darwin
}
set plat [regsub {Windows-amd64} $plat {Windows-x86}]
set plat [regsub {intel} $plat {x86}]
#package ifneeded tkpng 0.9 \
#    [list load [file join $dir tkpng-${plat}[info sharedlibextension]] tkpng]

package ifneeded tkpng 0.9 \
    [list loadTkPNGDll $dir $plat]
     
