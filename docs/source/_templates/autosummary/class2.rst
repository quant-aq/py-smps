{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

    {% block attributes %}
    {% if attributes %}

    ATTRIBUTES

    {% for item in attributes %}
    .. autoattribute:: {{ name }}.{{ item }}

    {% endfor %}
        
    {% endif %}
    {% endblock %}


    {% block methods %}
    {% if methods %}

    METHODS 

    {% for item in methods %}
    .. automethod:: {{ name }}.{{ item }}
    {% endfor %}
    {% endif %}
    {% endblock %}

    