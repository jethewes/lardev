art_make( MODULE_LIBRARIES nusimdata_SimulationBase
          ${LARRECO_LIB}
          ${LARDATA_LIB}
          ${ART_FRAMEWORK_CORE}
	  ${ART_FRAMEWORK_PRINCIPAL}
	  ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ${ART_ROOT_IO_TFILE_SUPPORT}
	  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
	  art_Persistency_Common canvas
	  art_Persistency_Provenance canvas
	  art_Utilities canvas
	  ${MF_MESSAGELOGGER}
	  ${MF_UTILITIES}
	  ${FHICLCPP}
	  cetlib cetlib_except
          ${ROOT_GEOM}
          ${ROOT_XMLIO}
          ${ROOT_GDML}
	  ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
install_scripts()
