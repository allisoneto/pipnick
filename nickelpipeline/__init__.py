
# # Adding a null handler to the package's root logger 
# # Prevents a default configuration being used if applications using this package don't set one:
# import logging
# logging.getLogger('nickelpipeline').addHandler(logging.NullHandler())