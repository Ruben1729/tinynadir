include_guard(GLOBAL)
# Directory groups organize source and documentation; selection stays per algorithm.
function(tinynadir_group_modules group result)
    if(group STREQUAL "core")
        set(_modules
            tinylinalg
            tinyrotation
            tinytime
            tinyeop
            tinyframes
            tinygeodesy
            tinyelements
            tinyanomaly
            tinystumpff
            tinyinterpolation
            tinypolynomial
            tinyelliptic
            tinyrk4
            tinyrkf45
        )
    elseif(group STREQUAL "orbit")
        set(_modules
            tinytle
            tinyomm
            tinysgp4
            tinykepler
            tinypkepler
            tinysecular
            tinymeanosculating
            tinycowell
        )
    elseif(group STREQUAL "environment")
        set(_modules
            tinysun
            tinymoon
            tinyephemeris
            tinyigrf
            tinygravity
            tinythirdbody
            tinyexponential
            tinymsis
            tinydtm
            tinyspaceweather
            tinydrag
            tinysrp
            tinyshadow
        )
    elseif(group STREQUAL "attitude")
        set(_modules
            tinycalibration
            tinytriad
            tinyquest
            tinydavenport
            tinygyro
            tinycomplementary
            tinymekf
        )
    elseif(group STREQUAL "control")
        set(_modules
            tinyguidance
            tinyslew
            tinybdot
            tinypd
            tinypid
            tinymagnetic
            tinywheels
            tinydesaturation
            tinyrigidbody
            tinydisturbance
        )
    elseif(group STREQUAL "mission")
        set(_modules
            tinyobservations
            tinygibbs
            tinyherrickgibbs
            tinygauss
            tinylaplace
            tinydoubler
            tinybatchls
            tinysgp4fit
            tinycovariance
            tinyunscented
            tinylambert
            tinyhohmann
            tinyplanechange
            tinyhills
            tinyrendezvous
            tinycoverage
        )
    elseif(group STREQUAL "all")
        set(_modules "")
        foreach(_group core orbit environment attitude control mission)
            tinynadir_group_modules(${_group} _members)
            list(APPEND _modules ${_members})
        endforeach()
    else()
        message(FATAL_ERROR "Unknown TinyNadir group: ${group}")
    endif()
    set(${result} "${_modules}" PARENT_SCOPE)
endfunction()

function(tinynadir_module_directory name result)
    foreach(_group core orbit environment attitude control mission)
        tinynadir_group_modules(${_group} _members)
        if(name IN_LIST _members)
            set(${result} "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../${_group}/${name}" PARENT_SCOPE)
            return()
        endif()
    endforeach()
    message(FATAL_ERROR "Unknown TinyNadir module: ${name}")
endfunction()
