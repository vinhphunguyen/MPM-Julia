module moduleMath
   type Vector2D   #node container
      fx::Real
      fy::Real

      function Vector2D()
         new(0.0, 0.0)
      end
      function Vector2D(f1::Real, f2::Real)
         new(f1, f2)
      end
   end

	type Vector3D   #node container
      f1::Real
      f2::Real
		f3::Real

      function Vector3D()
         new(0.0, 0.0, 0.0)
      end
      function Vector3D(f1::Real, f2::Real, f3::Real)
         new(f1, f2, f3)
      end
   end
end
